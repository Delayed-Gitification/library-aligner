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


def generate_illumina_se_data(tmp_path):
    random.seed(42)

    fasta_path = tmp_path / "reference.fasta"
    csv_path = tmp_path / "linkage.csv"
    fastq_path = tmp_path / "reads.fastq.gz"
    out_bam = tmp_path / "output.bam"

    # 1. Create the References (Standard genomic backbones)
    bases = ['A', 'C', 'G', 'T']
    ref_A = "".join(random.choices(bases, k=300))
    ref_B = "".join(random.choices(bases, k=300))

    with open(fasta_path, "w") as f:
        f.write(f">Genome_A\n{ref_A}\n")
        f.write(f">Genome_B\n{ref_B}\n")

    # 2. Assign barcodes in the CSV
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

    with gzip.open(fastq_path, "wt") as f:
        def write_read(qname, seq, expected_rname, strand):
            expected_qnames.append(qname)
            ground_truth[qname] = {"rname": expected_rname, "strand": strand}
            # Quality strings in Illumina are usually high (e.g., 'I' = Phred 40)
            f.write(f"@{qname}\n{seq}\n+\n{'I' * len(seq)}\n")

        configs = [
            (bc_A1, "Genome_A", ref_A),
            (bc_A2, "Genome_A", ref_A),
            (bc_B1, "Genome_B", ref_B),
            (bc_B2, "Genome_B", ref_B)
        ]

        read_counter = 0
        reads_per_bc = 100

        # Simulate a 5' inline barcode library prep:
        # 10bp Adapter + 15bp Barcode + 75bp Genomic Insert = 100bp Read
        adapter = "AATGATACGG"

        for bc, rname, ref_seq in configs:
            for _ in range(reads_per_bc // 2):
                read_counter += 1

                # Cut a random 75bp slice from the genome
                start_idx = random.randint(0, len(ref_seq) - 75)
                genomic_insert = ref_seq[start_idx: start_idx + 75]

                # Forward Read (Genomic insert is standard)
                read_fwd = adapter + bc + genomic_insert
                noisy_fwd = add_illumina_noise(read_fwd, error_rate=0.005)  # 1% noise
                write_read(f"read_{read_counter}_EXPECT_{rname}_STRAND_fwd", noisy_fwd, rname, "fwd")

                read_counter += 1

                # Reverse Read (Genomic insert is rev-comp'd, but 5' adapter & BC stay forward!)
                genomic_rev = revcomp(genomic_insert)
                read_rev = adapter + bc + genomic_rev
                noisy_rev = add_illumina_noise(read_rev, error_rate=0.005)
                write_read(f"read_{read_counter}_EXPECT_{rname}_STRAND_rev", noisy_rev, rname, "rev")

    return fasta_path, csv_path, fastq_path, out_bam, expected_qnames, ground_truth


# =============================================================================
# The Diagnostic Test
# =============================================================================

def test_illumina_single_end(tmp_path):
    fasta, csv, fastq, out_bam, expected_qnames, ground_truth = generate_illumina_se_data(tmp_path)

    import sys
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script_path = os.path.join(root_dir, "src", "library_aligner", "core.py")

    cmd = [
        sys.executable, script_path,
        "-fq", str(fastq),
        "-fa", str(fasta),
        "-c", str(csv),
        "-o", str(out_bam),
        "--barcode_start_pos", "10",  # <-- CRITICAL: Direct positional slicing
        "--barcode_end_pos", "25",  # <-- 15bp slice
        "--minimap2_preset", "splice:sr",  # <-- Requested preset
        "--minimum_distance", "1",  # Illumina is clean, 1 edit is plenty
        "--check_reverse_complement", "False"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        pytest.fail(f"Pipeline Crashed!\n{result.stderr}", pytrace=False)

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

            expected_rname = ground_truth[qname]["rname"]
            expected_strand = ground_truth[qname]["strand"]

            if read.reference_name == expected_rname:
                correct_templates += 1

            if expected_strand == "rev" and read.is_reverse:
                correct_strands += 1
            elif expected_strand == "fwd" and not read.is_reverse:
                correct_strands += 1

    missing_qnames = set(expected_qnames) - mapped_qnames
    mapped_val = len(mapped_qnames)
    target_val = len(expected_qnames) * 0.95

    # --- FAILURE REPORT ---
    if mapped_val < target_val or correct_templates != mapped_val or correct_strands != mapped_val:
        failed_reads_dump = ""
        dump_count = 0
        for name, seq, _ in mappy.fastx_read(str(fastq)):
            if name in missing_qnames and dump_count < 5:
                failed_reads_dump += f">{name}\n{seq}\n\n"
                dump_count += 1

        pipeline_log = [line for line in result.stderr.split('\n') if "Pipeline complete" in line]
        log_summary = pipeline_log[0] if pipeline_log else "No summary line found."

        pytest.fail(
            f"\n\n--- ILLUMINA SINGLE-END TEST FAILED ---\n"
            f"Expected 400 reads. Recovered: {mapped_val}/400\n"
            f"Template Mismatches: {mapped_val - correct_templates}\n"
            f"Strand Orientation Errors: {mapped_val - correct_strands}\n\n"
            f"Pipeline Summary:\n{log_summary}\n\n"
            f"--- SAMPLE OF FAILED READS ---\n"
            f"{failed_reads_dump}"
            f"-------------------------------------\n",
            pytrace=False
        )

    # --- SUCCESS REPORT ---
    print("\n\n" + "=" * 70)
    print("🟢 ILLUMINA SINGLE-END TEST PASSED 🟢")
    print("=" * 70)
    print(f"✅ Configuration: 100bp Reads | Positional Extraction [10:25] | Preset: sr:splice")
    print(f"✅ Total reads expected: {len(expected_qnames)}")
    print(f"✅ Total reads rescued: {mapped_val} / {len(expected_qnames)}")
    print(f"✅ Multi-Mapping accuracy: {correct_templates} / {mapped_val} mapped to correct genome")
    print(f"✅ Strand awareness: {correct_strands} / {mapped_val} aligned correctly")
    print("=" * 70 + "\n")