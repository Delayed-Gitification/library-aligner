import os
import subprocess
import random
import gzip
import pysam
import pytest
import mappy


# =============================================================================
# Helpers: Data Generation & Mutation
# =============================================================================

def revcomp(seq: str) -> str:
    trans = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(trans)[::-1]


def add_illumina_noise(seq: str, error_rate: float) -> str:
    # Illumina errors are overwhelmingly substitutions, very few indels.
    # This ensures our rigid positional extraction doesn't suffer frameshifts.
    bases = ['A', 'C', 'G', 'T']
    noisy_seq = []
    for base in seq:
        if random.random() < error_rate:
            noisy_seq.append(random.choice([b for b in bases if b != base]))
        else:
            noisy_seq.append(base)
    return "".join(noisy_seq)

def generate_illumina_pe_data(tmp_path):
    random.seed(42)

    fasta_path = tmp_path / "reference_pe.fasta"
    csv_path = tmp_path / "linkage_pe.csv"
    fastq_r1_path = tmp_path / "reads_R1.fastq.gz"
    fastq_r2_path = tmp_path / "reads_R2.fastq.gz"
    out_bam = tmp_path / "output_pe.bam"

    # 1. Create References — need to be long enough for PE150 inserts
    bases = ['A', 'C', 'G', 'T']
    ref_A = "".join(random.choices(bases, k=600))
    ref_B = "".join(random.choices(bases, k=600))

    with open(fasta_path, "w") as f:
        f.write(f">Genome_A\n{ref_A}\n")
        f.write(f">Genome_B\n{ref_B}\n")

    # 2. Assign barcodes
    bc_A1 = "GATTACAGATTACAG"
    bc_A2 = "CCATCCATCCATCCA"
    bc_B1 = "TTAGCTTAGCTTAGC"
    bc_B2 = "GGAAGGAAGGAAGGA"

    with open(csv_path, "w") as f:
        f.write("barcode,rname\n")
        f.write(f"{bc_A1},Genome_A\n")
        f.write(f"{bc_A2},Genome_A\n")
        f.write(f"{bc_B1},Genome_B\n")
        f.write(f"{bc_B2},Genome_B\n")

    expected_qnames = []
    ground_truth = {}

    # PE150 layout:
    # R1: 10bp Adapter + 15bp Barcode + 125bp Genomic = 150bp
    # R2: 150bp Genomic (reverse complement of insert end) = 150bp
    # Insert size: 250bp (so R1 and R2 overlap slightly, typical PE150)
    adapter = "AATGATACGGCGACC"  # 15bp to keep R1 round numbers — shift bc to [15:30]
    insert_size = 250
    read_len = 150
    genomic_in_r1 = read_len - len(adapter) - len(bc_A1)  # = 125bp (same for all bcs, all are 15bp)

    with gzip.open(fastq_r1_path, "wt") as f1, gzip.open(fastq_r2_path, "wt") as f2:

        def write_pair(qname, r1_seq, r2_seq, expected_rname):
            expected_qnames.append(qname)
            ground_truth[qname] = {"rname": expected_rname}
            f1.write(f"@{qname}\n{r1_seq}\n+\n{'I' * len(r1_seq)}\n")
            f2.write(f"@{qname}\n{r2_seq}\n+\n{'I' * len(r2_seq)}\n")

        configs = [
            (bc_A1, "Genome_A", ref_A),
            (bc_A2, "Genome_A", ref_A),
            (bc_B1, "Genome_B", ref_B),
            (bc_B2, "Genome_B", ref_B),
        ]

        read_counter = 0
        reads_per_bc = 100

        for bc, rname, ref_seq in configs:
            for _ in range(reads_per_bc):
                read_counter += 1

                # Pick a random insert of `insert_size` bp from the genome
                max_start = len(ref_seq) - insert_size
                start_idx = random.randint(0, max_start)
                insert = ref_seq[start_idx: start_idx + insert_size]

                # R1: adapter + barcode + first 125bp of insert (forward strand)
                r1_genomic = insert[:genomic_in_r1]
                r1 = adapter + bc + r1_genomic

                # R2: last 150bp of insert, reverse complemented
                # (sequencer reads the opposite strand from the other end)
                r2_genomic = insert[insert_size - read_len:]
                r2 = revcomp(r2_genomic)

                r1_noisy = add_illumina_noise(r1, error_rate=0.005)
                r2_noisy = add_illumina_noise(r2, error_rate=0.005)

                write_pair(f"read_{read_counter}_EXPECT_{rname}", r1_noisy, r2_noisy, rname)

    # Barcode lives at R1[15:30]
    bc_start = len(adapter)               # 15
    bc_end = len(adapter) + len(bc_A1)    # 30

    return fasta_path, csv_path, fastq_r1_path, fastq_r2_path, out_bam, expected_qnames, ground_truth, bc_start, bc_end


def test_illumina_paired_end(tmp_path):
    fasta, csv, fastq_r1, fastq_r2, out_bam, expected_qnames, ground_truth, bc_start, bc_end = generate_illumina_pe_data(tmp_path)

    import sys
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script_path = os.path.join(root_dir, "src", "library_aligner", "core.py")

    cmd = [
        sys.executable, script_path,
        "-fq", str(fastq_r1),
        "-fq2", str(fastq_r2),
        "-fa", str(fasta),
        "-c", str(csv),
        "-o", str(out_bam),
        "--barcode_start_pos", str(bc_start),
        "--barcode_end_pos", str(bc_end),
        "--barcode_read", "1",
        "--minimap2_preset", "sr",
        "--minimum_distance", "1",
        "--check_reverse_complement", "False",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        pytest.fail(f"Pipeline Crashed!\n{result.stderr}", pytrace=False)

    pysam.index(str(out_bam))

    mapped_qnames = set()
    correct_templates = 0
    r1_seen = set()
    r2_seen = set()

    with pysam.AlignmentFile(out_bam, "rb") as bam:
        for read in bam.fetch():
            if read.is_secondary or read.is_supplementary:
                continue

            qname = read.query_name

            if read.is_read1:
                r1_seen.add(qname)
            elif read.is_read2:
                r2_seen.add(qname)

            mapped_qnames.add(qname)

            expected_rname = ground_truth[qname]["rname"]
            if read.reference_name == expected_rname:
                correct_templates += 1

    # Pairs where both R1 and R2 are present
    complete_pairs = r1_seen & r2_seen

    missing_qnames = set(expected_qnames) - mapped_qnames
    mapped_val = len(mapped_qnames)
    target_val = len(expected_qnames) * 0.95

    if mapped_val < target_val or correct_templates < mapped_val * 2 * 0.95:
        failed_reads_dump = ""
        dump_count = 0
        for name, seq, _ in mappy.fastx_read(str(fastq_r1)):
            if name in missing_qnames and dump_count < 5:
                failed_reads_dump += f">{name} (R1)\n{seq}\n\n"
                dump_count += 1

        pipeline_log = [line for line in result.stderr.split('\n') if "Pipeline complete" in line]
        log_summary = pipeline_log[0] if pipeline_log else "No summary line found."

        pytest.fail(
            f"\n\n--- ILLUMINA PAIRED-END TEST FAILED ---\n"
            f"Expected {len(expected_qnames)} read pairs. Recovered: {mapped_val}/{len(expected_qnames)}\n"
            f"Complete pairs (R1+R2 both present): {len(complete_pairs)}/{len(expected_qnames)}\n"
            f"Template mismatches: {mapped_val * 2 - correct_templates}\n\n"
            f"Pipeline Summary:\n{log_summary}\n\n"
            f"--- SAMPLE OF FAILED R1 READS ---\n"
            f"{failed_reads_dump}"
            f"----------------------------------------\n",
            pytrace=False
        )

    print("\n\n" + "=" * 70)
    print("🟢 ILLUMINA PAIRED-END TEST PASSED 🟢")
    print("=" * 70)
    print(f"✅ Configuration: PE150 | Barcode on R1[{bc_start}:{bc_end}] | Preset: sr")
    print(f"✅ Total pairs expected:       {len(expected_qnames)}")
    print(f"✅ Total pairs recovered:      {mapped_val} / {len(expected_qnames)}")
    print(f"✅ Complete pairs (R1+R2):     {len(complete_pairs)} / {len(expected_qnames)}")
    print(f"✅ Template accuracy:          {correct_templates} / {mapped_val * 2} alignments")
    print("=" * 70 + "\n")