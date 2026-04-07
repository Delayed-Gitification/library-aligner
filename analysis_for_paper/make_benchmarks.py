import os, gzip, random, subprocess, time
import psutil, pysam, dnaio, pandas as pd

LA_PYTHON = "/Users/ogw/miniconda3/envs/library_aligner/bin/python"
LA_SCRIPT = "/Users/ogw/Documents/GitHub/library-aligner/src/library_aligner/core.py"
OUT_ROOT  = "/Users/ogw/Documents/GitHub/library-aligner/analysis_for_paper/benchmark_data"

WILDCARD = "NNNNNNNNNNNNNNN"   # 15 bp wildcard
F5       = "ACGTACGTAC" + "ACTGAC"        # 10 bp left flank
F3       = "TGCATGCATG" + "ATTAGG"        # 10 bp right flank
BASES    = ['A', 'C', 'G', 'T']

# Maximum reads sent to each aligner per condition.
# Keeping this low makes naive minimap2 tractable at large library sizes
# while still giving a fair comparison. Change here to affect all experiments.
MAX_READS = 1000

print("Config OK")


# ── Sequence utilities ────────────────────────────────────────────────────────

def mutate_seq(seq, error_rate, rng):
    """Introduce substitutions at a given rate."""
    return "".join(
        rng.choice([b for b in BASES if b != c]) if rng.random() < error_rate else c
        for c in seq
    )


# ── Benchmark data generators ─────────────────────────────────────────────────

def generate_benchmark_data(
    out_dir, num_barcodes, read_depth=100, seed=42,
    insert_error_rate=0.05, insert_length=400,
):
    """
    Generate three benchmark files for a synthetic barcode library:
      1. reference.fasta  — one entry per barcode template, with a wildcard
                            placeholder (NNNNNNNNNNNNNNN) in place of the barcode
      2. linkage.csv      — maps each barcode sequence to its reference entry name
      3. reads.fastq.gz   — noiseless reads with the actual barcode substituted in

    insert_error_rate controls per-base divergence between templates, mimicking
    the slight sequence variation across real library inserts.
    """
    rng = random.Random(seed)  # local RNG — doesn't affect global random state
    os.makedirs(out_dir, exist_ok=True)
    fasta_path = os.path.join(out_dir, "reference.fasta")
    csv_path   = os.path.join(out_dir, "linkage.csv")
    fastq_path = os.path.join(out_dir, "reads.fastq.gz")

    # Shared insert consensus; each template gets an independently mutated
    # copy, so insert sequences differ slightly across library members.
    insert_l = "".join(rng.choices(BASES, k=insert_length))
    insert_r = "".join(rng.choices(BASES, k=insert_length))
    barcodes = ["".join(rng.choices(BASES, k=15)) for _ in range(num_barcodes)]

    # Build wildcard templates: mut_insert_L + F5 + WILDCARD + F3 + mut_insert_R
    # The WILDCARD placeholder is what library-aligner indexes and extracts from.
    bc_templates = {}
    for i, bc in enumerate(barcodes):
        mut_l    = mutate_seq(insert_l, insert_error_rate, rng)
        mut_r    = mutate_seq(insert_r, insert_error_rate, rng)
        template = mut_l + F5 + WILDCARD + F3 + mut_r
        bc_templates[bc] = (f"Ref_{i}", template)

    # FASTA holds wildcard templates (not real barcodes) — this is what
    # library-aligner indexes. CSV links each barcode to its reference entry.
    with open(fasta_path, "w") as fa, open(csv_path, "w") as csv_f:
        csv_f.write("barcode,rname\n")
        for bc, (rname, template) in bc_templates.items():
            fa.write(f">{rname}\n{template}\n")
            csv_f.write(f"{bc},{rname}\n")

    # Emit read_depth perfect reads per barcode (WILDCARD replaced with actual
    # barcode). ground_truth is the per-read oracle used later by
    # measure_accuracy_la / measure_accuracy_naive.
    ground_truth = {}
    with gzip.open(fastq_path, "wt") as fq:
        for i, (bc, (rname, template)) in enumerate(bc_templates.items()):
            resolved = template.replace(WILDCARD, bc)
            for j in range(read_depth):
                qname = f"read_{i}_{j}"
                ground_truth[qname] = bc
                fq.write(f"@{qname}\n{resolved}\n+\n{'I' * len(resolved)}\n")

    return fasta_path, csv_path, fastq_path, barcodes, ground_truth


def generate_naive_mrna_fasta(out_dir, fastq_path, barcodes):
    """
    Build a FASTA of spliced mRNA sequences (entry names = barcodes) for use
    as a badread simulation source in spliced experiments. Unlike naive_fasta
    which is built from the genomic reference, this is built directly from the
    noiseless mRNA reads written by generate_spliced_benchmark.
    """
    mrna_fasta_path = os.path.join(out_dir, "naive_mrna_reference.fasta")
    bc_set = set(barcodes)
    seen   = {}
    with dnaio.open(fastq_path) as f:
        for record in f:
            qname = record.name.split()[0]
            # read names are read_{i}_{j} — take one representative per barcode
            parts = qname.split("_")
            bc_idx = int(parts[1])
            bc = barcodes[bc_idx]
            if bc not in seen:
                seen[bc] = record.sequence
    with open(mrna_fasta_path, "w") as fa:
        for bc, seq in seen.items():
            fa.write(f">{bc}\n{seq}\n")
    return mrna_fasta_path


def generate_naive_fasta(out_dir, fasta_path, csv_path):
    """
    Build a concatenated FASTA where every barcode gets its own entry
    (wildcard replaced with the actual barcode). Entry names ARE barcodes —
    this is both the naive alignment reference and the simulation source.
    """
    naive_fasta_path = os.path.join(out_dir, "naive_reference.fasta")
    references, current_name = {}, None
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_name = line[1:]
                references[current_name] = ""
            else:
                references[current_name] += line

    df = pd.read_csv(csv_path)
    with open(naive_fasta_path, "w") as fa:
        for _, row in df.iterrows():
            bc       = row["barcode"]
            template = references[row["rname"]]
            fa.write(f">{bc}\n{template.replace(WILDCARD, bc)}\n")
    return naive_fasta_path


def generate_spliced_benchmark(
    out_dir, num_barcodes, read_depth=100, seed=42,
    insert_error_rate=0.05,
    exon_lengths=(200, 200, 200),
    intron_length=200,
):
    """
    Generate benchmark data for spliced alignment.
    Reference FASTA = pre-mRNA (exons + canonical GT-AG introns).
    Reads          = spliced mRNA (exons only, barcode in exon 1).
    This tests that the splice preset correctly handles reads that span
    intron-exon boundaries.
    """
    rng = random.Random(seed)
    os.makedirs(out_dir, exist_ok=True)
    fasta_path = os.path.join(out_dir, "reference.fasta")
    csv_path   = os.path.join(out_dir, "linkage.csv")
    fastq_path = os.path.join(out_dir, "reads.fastq.gz")

    barcodes     = ["".join(rng.choices(BASES, k=15)) for _ in range(num_barcodes)]
    n_exons      = len(exon_lengths)
    ground_truth = {}

    with open(fasta_path, "w") as fa, open(csv_path, "w") as csv_f, \
         gzip.open(fastq_path, "wt") as fq:
        csv_f.write("barcode,rname\n")

        for i, bc in enumerate(barcodes):
            rname = f"Ref_{i}"

            # Each exon gets independently mutated sequence, simulating
            # per-member insert divergence across the library.
            exons = [
                mutate_seq("".join(rng.choices(BASES, k=l)), insert_error_rate, rng)
                for l in exon_lengths
            ]

            # Canonical GT-AG splice sites with random intronic sequence between them
            introns = [
                "GT" + "".join(rng.choices(BASES, k=intron_length - 4)) + "AG"
                for _ in range(n_exons - 1)
            ]

            # Barcode placed in exon 1: left + F5 + WILDCARD/bc + F3 + right
            # Split at 20 bp so the barcode sits within the exon body, not at the very edge
            split       = min(20, len(exons[0]))
            left, right = exons[0][:split], exons[0][split:]
            exon1_wc    = left + F5 + WILDCARD + F3 + right
            exon1_bc    = left + F5 + bc       + F3 + right

            # Genomic reference written to FASTA: exons interleaved with introns
            genomic = exon1_wc + "".join(
                intron + exon for intron, exon in zip(introns, exons[1:])
            )

            # mRNA read: introns removed, actual barcode substituted.
            # The aligner must correctly skip the introns to assign the right barcode.
            mrna = exon1_bc + "".join(exons[1:])

            fa.write(f">{rname}\n{genomic}\n")
            csv_f.write(f"{bc},{rname}\n")

            for j in range(read_depth):
                qname = f"read_{i}_{j}"
                ground_truth[qname] = bc
                fq.write(f"@{qname}\n{mrna}\n+\n{'I' * len(mrna)}\n")

    return fasta_path, csv_path, fastq_path, barcodes, ground_truth


# ── Subsampling ───────────────────────────────────────────────────────────────

def subsample_fastq(fastq_path, ground_truth, max_reads, out_dir, seed=42):
    """
    Subsample a FASTQ to at most max_reads reads, guided by the ground truth.
    Returns (new_fastq_path, subsampled_ground_truth).
    If len(ground_truth) <= max_reads the inputs are returned unchanged.
    """
    rng        = random.Random(seed)
    all_qnames = list(ground_truth.keys())

    if len(all_qnames) <= max_reads:
        return fastq_path, ground_truth

    kept     = set(rng.sample(all_qnames, max_reads))
    sub_gt   = {q: ground_truth[q] for q in kept}
    out_path = os.path.join(out_dir, "subsampled.fastq.gz")

    with dnaio.open(fastq_path) as f_in, gzip.open(out_path, "wt") as f_out:
        for record in f_in:
            qname = record.name.split()[0]
            if qname in kept:
                f_out.write(
                    f"@{qname}\n{record.sequence}\n+\n{record.qualities}\n"
                )

    return out_path, sub_gt


# ── Noisy read simulators ─────────────────────────────────────────────────────

def simulate_badread(naive_fasta_path, out_dir, quantity="5x", seed=42):
    """
    Simulate ONT reads with Badread (https://github.com/rrwick/Badread).
    Uses the naive FASTA (entry names = barcodes) so ground truth is free.
    Requires: badread on PATH.
    """
    fastq_path = os.path.join(out_dir, "badread_reads.fastq.gz")
    cmd = [
        "badread", "simulate",
        "--reference",    naive_fasta_path,
        "--quantity",     quantity,
        "--error_model",  "nanopore2023",
        "--qscore_model", "nanopore2023",
        "--seed",         str(seed),
        "--length",       "1000,100",
    ]
    proc = subprocess.run(cmd, capture_output=True)
    with gzip.open(fastq_path, "wt") as out:
        out.write(proc.stdout.decode())
    return fastq_path


def simulate_art(naive_fasta_path, out_dir, read_length=150, fold_coverage=10, seed=42):
    """
    Simulate Illumina reads with ART (https://www.niehs.nih.gov/research/resources/software/biostatistics/art).
    Uses the naive FASTA (entry names = barcodes).
    Requires: art_illumina on PATH.
    """
    prefix = os.path.join(out_dir, "art_reads")
    cmd = [
        "art_illumina",
        "-ss", "HS25",
        "-i",  naive_fasta_path,
        "-l",  str(read_length),
        "-f",  str(fold_coverage),
        "-o",  prefix,
        "-na",           # no alignment output
        "-rs", str(seed),
        "-q",            # quiet
    ]
    subprocess.run(cmd, capture_output=True)
    art_fq     = prefix + ".fq"
    fastq_path = prefix + ".fastq.gz"
    with open(art_fq) as f_in, gzip.open(fastq_path, "wt") as f_out:
        f_out.write(f_in.read())
    os.remove(art_fq)
    return fastq_path


def parse_gt_badread(fastq_path):
    """
    Parse ground truth from Badread read names.
    Format: '<uuid> <entry_name>,<start>-<end>,<strand> ...'
    Entry name in our naive FASTA IS the barcode, so gt is free.
    """
    gt = {}
    with dnaio.open(fastq_path) as f:
        for record in f:
            qname = record.name.split()[0]
            parts = record.name.split()
            if len(parts) > 1:
                bc = parts[1].split(",")[0]
                gt[qname] = bc
    return gt


def parse_gt_art(fastq_path, known_barcodes):
    """
    Parse ground truth from ART read names.
    Format: '<entry_name>-<readnumber>' where entry_name IS the barcode.
    Sort barcodes by length descending to avoid prefix collisions.
    """
    sorted_bcs = sorted(known_barcodes, key=len, reverse=True)
    gt = {}
    with dnaio.open(fastq_path) as f:
        for record in f:
            qname = record.name.split()[0]
            for bc in sorted_bcs:
                if qname.startswith(bc + "-") or qname.startswith(bc + "/"):
                    gt[qname] = bc
                    break
    return gt


# ── Process monitoring ────────────────────────────────────────────────────────

def _poll_process(proc):
    """Poll a subprocess, returning (elapsed_s, peak_rss_MB)."""
    ps_proc  = psutil.Process(proc.pid)
    peak_mem = 0
    start    = time.time()
    while proc.poll() is None:
        try:
            mem = ps_proc.memory_info().rss
            if mem > peak_mem:
                peak_mem = mem
        except psutil.NoSuchProcess:
            break
        time.sleep(0.05)
    return time.time() - start, peak_mem / 1024**2


# ── Alignment benchmarks ──────────────────────────────────────────────────────

def benchmark_library_aligner(
    fasta_path, csv_path, fastq_path, out_dir,
    preset="splice", kmer=10, w=4, no_cache=False, barcode_max_edits=0,
):
    tag     = "nocache" if no_cache else "cached"
    out_bam = os.path.join(out_dir, f"la_{tag}.bam")
    cmd = [
        LA_PYTHON, LA_SCRIPT,
        "-fq", fastq_path, "-fa", fasta_path, "-c",  csv_path,
        "-o",  out_bam,    "-wc", WILDCARD,
        "--minimap2_preset", preset,
        "--minimap2_kmer",   str(kmer),
        "--minimap2_w",      str(w),
        "--max_bc_edit_distance", str(barcode_max_edits),
    ]
    if no_cache:
        cmd.append("--no_mappy_cache")
    elapsed, peak_mem = _poll_process(subprocess.Popen(cmd))
    return elapsed, peak_mem, out_bam


def benchmark_naive(naive_fasta_path, fastq_path, out_dir, preset="splice", kmer=10, w=4):
    out_sam = os.path.join(out_dir, "naive_output.sam")
    cmd = [
        "minimap2", "-ax", preset,
        "-k", str(kmer), "-w", str(w),
        naive_fasta_path, fastq_path,
        "-o", out_sam,
    ]
    elapsed, peak_mem = _poll_process(subprocess.Popen(cmd))
    return elapsed, peak_mem, out_sam


# ── Accuracy measurement ──────────────────────────────────────────────────────

def measure_accuracy_la(bam_path, ground_truth):
    """Check BC:Z tag in library-aligner output against ground truth."""
    total = correct = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch():
            if read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            if qname not in ground_truth:
                continue
            total += 1
            if read.has_tag("BC") and read.get_tag("BC") == ground_truth[qname]:
                correct += 1
    return correct / total if total > 0 else 0.0


def measure_accuracy_naive(sam_path, ground_truth):
    """
    Naive FASTA entry names ARE barcodes, so check whether the reference
    the read mapped to matches the expected barcode in ground_truth.
    """
    total = correct = 0
    with pysam.AlignmentFile(sam_path, "rb") as sam:
        for read in sam.fetch():
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            qname = read.query_name
            if qname not in ground_truth:
                continue
            total += 1
            if read.reference_name == ground_truth[qname]:
                correct += 1
    return correct / total if total > 0 else 0.0


# ── Experiment runner ─────────────────────────────────────────────────────────

def run_experiment(
    label, scales, fixed_kwargs, vary_key, out_root,
    preset="splice", kmer=10, w=4,
    data_fn=generate_benchmark_data,
    simulator=None, barcode_max_edits=0,
):
    """
    Run a benchmark experiment across multiple values of a single parameter
    (vary_key), keeping all other parameters fixed. For example, vary_key="num_barcodes"
    with scales=[10, 100, 1000] will run three conditions testing how performance
    changes with library size.

    data_fn controls what kind of synthetic data is generated — pass
    generate_spliced_benchmark to test spliced alignment instead of the default.

    simulator optionally layers realistic read errors on top of the generated
    references: "badread" for ONT noise, "art" for Illumina noise. When set,
    noisy reads are simulated from the naive FASTA and ground truth is parsed
    directly from the simulator's read names (no manual labelling needed, since
    naive FASTA entry names ARE barcodes).

    All three aligners (library-aligner cached, no-cache, naive minimap2) are
    run on the same subsampled reads for each condition. Saves results.csv under
    out_root/label/ and returns the dataframe.
    """
    rows = []
    for val in scales:
        kwargs  = {**fixed_kwargs, vary_key: val}
        out_dir = os.path.join(out_root, label, str(val))
        print(f"  [{label}] {vary_key}={val} ...")

        fasta, csv, fastq, barcodes, gt = data_fn(out_dir, **kwargs)
        naive_fasta = generate_naive_fasta(out_dir, fasta, csv)

        if simulator is not None:
            # Replace the noiseless fastq with simulator output.
            if simulator == "badread":
                target_reads = MAX_READS * 2
                sim_fasta = generate_naive_mrna_fasta(out_dir, fastq, barcodes) if data_fn == generate_spliced_benchmark else naive_fasta
                with open(sim_fasta) as f:
                    actual_bases = sum(len(l.strip()) for l in f if not l.startswith(">"))
                quantity = f"{max(1, round((target_reads * 1000) / actual_bases))}x"
                fastq = simulate_badread(sim_fasta, out_dir, quantity=quantity)
                gt    = parse_gt_badread(fastq)
            elif simulator == "art":
                fastq = simulate_art(naive_fasta, out_dir)
                gt    = parse_gt_art(fastq, barcodes)
            else:
                raise ValueError(f"Unknown simulator: {simulator}")
            print(f"    {simulator} generated {len(gt)} reads")

        # Cap reads sent to each aligner so naive minimap2 stays tractable at
        # large library sizes while keeping the comparison fair between methods.
        fastq_sub, gt_sub = subsample_fastq(fastq, gt, MAX_READS, out_dir)

        # Run all three aligner configurations on the same subsampled reads
        la_t,  la_m,  la_bam  = benchmark_library_aligner(fasta, csv, fastq_sub, out_dir, preset, kmer, w, barcode_max_edits=barcode_max_edits)
        lac_t, lac_m, lac_bam = benchmark_library_aligner(fasta, csv, fastq_sub, out_dir, preset, kmer, w, no_cache=True, barcode_max_edits=barcode_max_edits)
        nv_t,  nv_m,  nv_sam  = benchmark_naive(naive_fasta, fastq_sub, out_dir, preset, kmer, w)

        la_acc  = measure_accuracy_la(la_bam,  gt_sub)
        lac_acc = measure_accuracy_la(lac_bam, gt_sub)
        nv_acc  = measure_accuracy_naive(nv_sam, gt_sub)

        rows.append({
            vary_key:              val,
            "la_time":             la_t,  "la_mem":         la_m,  "la_accuracy":         la_acc,
            "la_nocache_time":     lac_t, "la_nocache_mem": lac_m, "la_nocache_accuracy": lac_acc,
            "naive_time":          nv_t,  "naive_mem":      nv_m,  "naive_accuracy":      nv_acc,
        })

    df = pd.DataFrame(rows)
    os.makedirs(os.path.join(out_root, label), exist_ok=True)
    df.to_csv(os.path.join(out_root, label, "results.csv"), index=False)
    return df


print("Helper functions defined.")


# ── Experiments ───────────────────────────────────────────────────────────────

BARCODE_SCALES = [500]
INSERT_LEN     = 400
ERROR_RATE     = 0.02

# # Exp 1: how does performance scale with library size?
# df_exp1 = run_experiment(
#     label        = "exp1",
#     scales       = BARCODE_SCALES,
#     fixed_kwargs = dict(insert_length=INSERT_LEN, insert_error_rate=ERROR_RATE, read_depth=200),
#     vary_key     = "num_barcodes",
#     out_root     = OUT_ROOT,
# )
# print(df_exp1.to_string(index=False))

# # Exp 2: how does performance scale with insert length?
# df_exp2 = run_experiment(
#     label        = "exp2",
#     scales       = [50, 200, 400, 1000, 2000, 3000],
#     fixed_kwargs = dict(num_barcodes=250, read_depth=20, insert_error_rate=ERROR_RATE),
#     vary_key     = "insert_length",
#     out_root     = OUT_ROOT,
# )
# print(df_exp2.to_string(index=False))

# # Exp 3: how does performance degrade as inserts diverge from each other?
# df_exp3 = run_experiment(
#     label        = "exp3",
#     scales       = [0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0],
#     fixed_kwargs = dict(num_barcodes=500, read_depth=4, insert_length=INSERT_LEN),
#     vary_key     = "insert_error_rate",
#     out_root     = OUT_ROOT,
# )
# print(df_exp3.to_string(index=False))

# # Exp 4: library size scaling with realistic ONT noise (Badread)
# df_exp4 = run_experiment(
#     label        = "exp4_ont",
#     scales       = BARCODE_SCALES,
#     fixed_kwargs = dict(insert_length=INSERT_LEN, insert_error_rate=ERROR_RATE, read_depth=20),
#     vary_key     = "num_barcodes",
#     out_root     = OUT_ROOT,
#     simulator    = "badread",
#     preset       = "splice",
#     kmer         = 10,
#     w            = 4,
# )
# print(df_exp4.to_string(index=False))

# # Exp 5: library size scaling with realistic Illumina noise (ART).
# # Short insert_length=100 keeps reads within a single fragment for short-read alignment.
# df_exp5 = run_experiment(
#     label        = "exp5_illumina",
#     scales       = BARCODE_SCALES,
#     fixed_kwargs = dict(insert_length=100, insert_error_rate=ERROR_RATE, read_depth=20),
#     vary_key     = "num_barcodes",
#     out_root     = OUT_ROOT,
#     simulator    = "art",
#     preset       = "splice:sr",
#     kmer         = 10,
#     w            = 4,
# )
# print(df_exp5.to_string(index=False))

# # Exp 6: library size scaling with spliced reads aligned to genomic references
# df_exp6 = run_experiment(
#     label        = "exp6_spliced",
#     scales       = BARCODE_SCALES,
#     fixed_kwargs = dict(insert_error_rate=ERROR_RATE, read_depth=20,
#                         exon_lengths=(200, 200, 200), intron_length=200),
#     vary_key     = "num_barcodes",
#     out_root     = OUT_ROOT,
#     data_fn      = generate_spliced_benchmark,
#     preset       = "splice",
#     kmer         = 10,
#     w            = 4,
# )
# print(df_exp6.to_string(index=False))

# Exp 7: worst-case — spliced reads, short inserts, ONT noise
df_exp7 = run_experiment(
    label        = "exp7_spliced_ont",
    scales       = BARCODE_SCALES,
    fixed_kwargs = dict(insert_error_rate=ERROR_RATE, read_depth=20,
                        exon_lengths=(40, 50), intron_length=200),
    vary_key     = "num_barcodes",
    out_root     = OUT_ROOT,
    data_fn      = generate_spliced_benchmark,
    simulator    = "badread",
    preset       = "splice",
    kmer         = 7,
    w            = 3,
    barcode_max_edits = 5,
)
print(df_exp7.to_string(index=False))