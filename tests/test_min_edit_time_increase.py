import os
import subprocess
import random
import gzip
import time
import pysam
import pytest
import mappy


# =============================================================================
# Helpers
# =============================================================================

def revcomp(seq: str) -> str:
    trans = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(trans)[::-1]


def add_sequencing_noise(seq: str, error_rate: float) -> str:
    bases = ['A', 'C', 'G', 'T']
    return "".join(
        random.choice([b for b in bases if b != base]) if random.random() < error_rate else base
        for base in seq
    )


def generate_distinct_barcodes(n: int, length: int = 15, min_edit_dist: int = 4) -> list[str]:
    """
    Generate N barcodes that are all at least min_edit_dist apart from each other.
    This ensures fuzzy matching at distance 1 can never be ambiguous.
    """
    import edlib
    bases = ['A', 'C', 'G', 'T']
    accepted = []
    attempts = 0

    while len(accepted) < n:
        attempts += 1
        if attempts > 100_000:
            raise RuntimeError(
                f"Could not generate {n} barcodes with min edit distance {min_edit_dist} "
                f"after 100k attempts. Try fewer barcodes or a shorter min_edit_dist."
            )
        candidate = "".join(random.choices(bases, k=length))
        if all(
            edlib.align(candidate, bc, mode="NW", task="distance")["editDistance"] >= min_edit_dist
            for bc in accepted
        ):
            accepted.append(candidate)

    return accepted


def generate_large_library_data(tmp_path, n_plasmids: int = 2000, reads_per_barcode: int = 50):
    random.seed(99)

    fasta_path = tmp_path / "reference_large.fasta"
    csv_path = tmp_path / "linkage_large.csv"
    fastq_path = tmp_path / "reads_large.fastq.gz"
    out_bam_strict = tmp_path / "output_strict.bam"   # min_dist=0
    out_bam_fuzzy  = tmp_path / "output_fuzzy.bam"    # min_dist=1

    f5 = "ACGTACGTAC"   # 10bp
    f3 = "TGCATGCATG"   # 10bp
    wildcard = "NNNNNNNNNNNNNNN"  # 15bp

    # Each plasmid has a unique random backbone so alignments are unambiguous
    bases = ['A', 'C', 'G', 'T']

    print(f"Generating {n_plasmids} distinct barcodes (min pairwise edit dist = 4)...")
    barcodes = generate_distinct_barcodes(n_plasmids, length=15, min_edit_dist=4)

    # Write FASTA — one template per plasmid, wildcard in the middle
    backbone_len = 200
    plasmid_seqs = {}
    with open(fasta_path, "w") as f:
        for i, bc in enumerate(barcodes):
            rname = f"Plasmid_{i+1}"
            left  = "".join(random.choices(bases, k=backbone_len))
            right = "".join(random.choices(bases, k=backbone_len))
            template = left + f5 + wildcard + f3 + right
            plasmid_seqs[bc] = (rname, template)
            f.write(f">{rname}\n{template}\n")

    # CSV linkage
    with open(csv_path, "w") as f:
        f.write("barcode,rname\n")
        for bc, (rname, _) in plasmid_seqs.items():
            f.write(f"{bc},{rname}\n")

    expected_qnames = []
    ground_truth = {}   # qname -> {"rname": ..., "strand": ...}

    total_reads = n_plasmids * reads_per_barcode
    print(f"Writing {total_reads:,} reads ({reads_per_barcode} per barcode × {n_plasmids} barcodes)...")

    with gzip.open(fastq_path, "wt") as fq:
        def write_read(qname, seq, rname, strand):
            expected_qnames.append(qname)
            ground_truth[qname] = {"rname": rname, "strand": strand}
            fq.write(f"@{qname}\n{seq}\n+\n{'I' * len(seq)}\n")

        read_counter = 0
        for bc, (rname, template) in plasmid_seqs.items():
            resolved = template.replace(wildcard, bc)
            bc_idx = resolved.find(bc)
            read_start = bc_idx - 50  # centre the barcode in the read

            for _ in range(reads_per_barcode // 2):
                read_counter += 1
                raw = resolved[read_start: read_start + 150]

                # Forward — 1% noise, barcode region may get 1 error occasionally
                fwd = add_sequencing_noise(raw, error_rate=0.01)
                write_read(f"read_{read_counter}_fwd_{rname}", fwd, rname, "fwd")

                read_counter += 1

                # Reverse — same noise level
                rev = add_sequencing_noise(revcomp(raw), error_rate=0.01)
                write_read(f"read_{read_counter}_rev_{rname}", rev, rname, "rev")

    return fasta_path, csv_path, fastq_path, out_bam_strict, out_bam_fuzzy, expected_qnames, ground_truth


# =============================================================================
# Shared runner
# =============================================================================

def run_pipeline(script_path, fastq, fasta, csv, out_bam, min_dist: int) -> tuple[subprocess.CompletedProcess, float]:
    import sys
    cmd = [
        sys.executable, script_path,
        "-fq", str(fastq),
        "-fa", str(fasta),
        "-c", str(csv),
        "-o", str(out_bam),
        "-wc", "NNNNNNNNNNNNNNN",
        "--five_prime_flank_length", "10",
        "--three_prime_flank_length", "10",
        "--minimum_distance", str(min_dist),
        "--minimap2_preset", "splice",
    ]
    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    return result, elapsed


def score_bam(out_bam, ground_truth) -> dict:
    pysam.index(str(out_bam))
    mapped_qnames = set()
    correct_templates = 0
    correct_strands = 0

    with pysam.AlignmentFile(out_bam, "rb") as bam:
        for read in bam.fetch():
            if read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            mapped_qnames.add(qname)
            gt = ground_truth[qname]
            if read.reference_name == gt["rname"]:
                correct_templates += 1
            if gt["strand"] == "rev" and read.is_reverse:
                correct_strands += 1
            elif gt["strand"] == "fwd" and not read.is_reverse:
                correct_strands += 1

    n = len(mapped_qnames)
    return {
        "mapped": n,
        "correct_templates": correct_templates,
        "correct_strands": correct_strands,
        "template_accuracy": correct_templates / n if n else 0,
        "strand_accuracy": correct_strands / n if n else 0,
    }


# =============================================================================
# The Test
# =============================================================================

N_PLASMIDS = 200
READS_PER_BARCODE = 50   # 200 × 50 = 10,000 reads total


def test_large_library_strict_vs_fuzzy(tmp_path):
    fasta, csv, fastq, out_bam_strict, out_bam_fuzzy, expected_qnames, ground_truth = \
        generate_large_library_data(tmp_path, n_plasmids=N_PLASMIDS, reads_per_barcode=READS_PER_BARCODE)

    import sys
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script_path = os.path.join(root_dir, "src", "library_aligner", "core.py")

    total_reads = len(expected_qnames)

    # --- Run strict (min_dist=0) ---
    print(f"\nRunning strict pass (min_dist=0) over {total_reads:,} reads...")
    result_strict, t_strict = run_pipeline(script_path, fastq, fasta, csv, out_bam_strict, min_dist=0)
    if result_strict.returncode != 0:
        pytest.fail(f"Strict pipeline crashed!\n{result_strict.stderr}", pytrace=False)

    # --- Run fuzzy (min_dist=1) ---
    print(f"Running fuzzy pass (min_dist=1) over {total_reads:,} reads...")
    result_fuzzy, t_fuzzy = run_pipeline(script_path, fastq, fasta, csv, out_bam_fuzzy, min_dist=2)
    if result_fuzzy.returncode != 0:
        pytest.fail(f"Fuzzy pipeline crashed!\n{result_fuzzy.stderr}", pytrace=False)

    # --- Score both ---
    strict = score_bam(out_bam_strict, ground_truth)
    fuzzy  = score_bam(out_bam_fuzzy,  ground_truth)

    recovered_by_fuzzy = fuzzy["mapped"] - strict["mapped"]
    miscalls = fuzzy["mapped"] - fuzzy["correct_templates"]

    # --- Correctness guard: fuzzy must not introduce template miscalls ---
    if fuzzy["template_accuracy"] < 0.999:
        pytest.fail(
            f"Fuzzy matching introduced template miscalls!\n"
            f"Template accuracy: {fuzzy['correct_templates']}/{fuzzy['mapped']} "
            f"({fuzzy['template_accuracy']:.3%})",
            pytrace=False
        )

    # --- Print comparison report ---
    print("\n\n" + "=" * 70)
    print("🟢 LARGE LIBRARY STRICT vs FUZZY COMPARISON PASSED 🟢")
    print("=" * 70)
    print(f"  Library size:        {N_PLASMIDS} plasmids × {READS_PER_BARCODE} reads = {total_reads:,} total")
    print(f"  Noise level:         1% substitution rate")
    print(f"  Barcode separation:  ≥4 edits between any two barcodes")
    print()
    print(f"  {'Metric':<30} {'Strict (dist=0)':>18} {'Fuzzy (dist=1)':>18}")
    print(f"  {'-'*66}")
    print(f"  {'Reads mapped':<30} {strict['mapped']:>17,} {fuzzy['mapped']:>17,}")
    print(f"  {'Recovery rate':<30} {strict['mapped']/total_reads:>17.1%} {fuzzy['mapped']/total_reads:>17.1%}")
    print(f"  {'Template accuracy':<30} {strict['template_accuracy']:>17.3%} {fuzzy['template_accuracy']:>17.3%}")
    print(f"  {'Strand accuracy':<30} {strict['strand_accuracy']:>17.3%} {fuzzy['strand_accuracy']:>17.3%}")
    print(f"  {'Runtime (s)':<30} {t_strict:>17.1f} {t_fuzzy:>17.1f}")
    print(f"  {'-'*66}")
    print(f"  Reads rescued by fuzzy:      {recovered_by_fuzzy:,}")
    print(f"  Miscalls introduced:         {miscalls}")
    print(f"  Runtime overhead:            +{t_fuzzy - t_strict:.1f}s ({(t_fuzzy/t_strict - 1)*100:.0f}%)")
    print("=" * 70 + "\n")