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


def add_sequencing_noise(seq: str, error_rate: float) -> str:
    bases = ['A', 'C', 'G', 'T']
    noisy_seq = []
    for base in seq:
        if random.random() < error_rate:
            noisy_seq.append(random.choice([b for b in bases if b != base]))
        else:
            noisy_seq.append(base)
    return "".join(noisy_seq)


def generate_slightly_noisy_data(tmp_path):
    """Generates 4 reads with highly visible sequences, but injects 2% random noise."""
    random.seed(42)  # Ensures the "random" noise is identical every time you run it

    fasta_path = tmp_path / "reference.fasta"
    csv_path = tmp_path / "linkage.csv"
    fastq_path = tmp_path / "reads.fastq.gz"
    out_bam = tmp_path / "output.bam"

    f5 = "ACGTACGTAC"  # 10bp
    f3 = "TGCATGCATG"  # 10bp
    wildcard = "NNNNNNNNNNNNNNN"

    # The backbone is pure As and Ts so the barcode pops out visually
    plasmid_seq = ("A" * 100) + f5 + wildcard + f3 + ("T" * 100)

    with open(fasta_path, "w") as f:
        f.write(f">Plasmid_1\n{plasmid_seq}\n")

    # We will use two highly distinct barcodes
    bc_1 = "GATTACAGATTACAG"
    bc_2 = "CCATCCATCCATCCA"

    with open(csv_path, "w") as f:
        f.write("barcode,rname\n")
        f.write(f"{bc_1},Plasmid_1\n")
        f.write(f"{bc_2},Plasmid_1\n")

    expected_qnames = []

    with gzip.open(fastq_path, "wt") as f:
        def write_read(qname, seq):
            expected_qnames.append(qname)
            f.write(f"@{qname}\n{seq}\n+\n{'I' * len(seq)}\n")

        # Create the true template for Barcode 1
        template_1 = plasmid_seq.replace(wildcard, bc_1)
        idx_1 = template_1.find(bc_1)

        # Read 1: Noisy Forward BC1 (150bp)
        seq_fwd_1 = template_1[idx_1 - 50: idx_1 + 100]
        noisy_fwd_1 = add_sequencing_noise(seq_fwd_1, error_rate=0.02)
        write_read("read_1_fwd_bc1_noisy", noisy_fwd_1)

        # Read 2: Noisy Reverse BC1 (150bp)
        seq_rev_1 = revcomp(template_1[idx_1 - 50: idx_1 + 100])
        noisy_rev_1 = add_sequencing_noise(seq_rev_1, error_rate=0.02)
        write_read("read_2_rev_bc1_noisy", noisy_rev_1)

        # Create the true template for Barcode 2
        template_2 = plasmid_seq.replace(wildcard, bc_2)
        idx_2 = template_2.find(bc_2)

        # Read 3: Noisy Forward BC2 (150bp)
        seq_fwd_2 = template_2[idx_2 - 50: idx_2 + 100]
        noisy_fwd_2 = add_sequencing_noise(seq_fwd_2, error_rate=0.02)
        write_read("read_3_fwd_bc2_noisy", noisy_fwd_2)

        # Read 4: Noisy Reverse BC2 (150bp)
        seq_rev_2 = revcomp(template_2[idx_2 - 50: idx_2 + 100])
        noisy_rev_2 = add_sequencing_noise(seq_rev_2, error_rate=0.02)
        write_read("read_4_rev_bc2_noisy", noisy_rev_2)

    return fasta_path, csv_path, fastq_path, out_bam, expected_qnames


# =============================================================================
# The Diagnostic Test
# =============================================================================

def test_wildcard_slight_noise(tmp_path):
    fasta, csv, fastq, out_bam, expected_qnames = generate_slightly_noisy_data(tmp_path)

    import sys
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script_path = os.path.join(root_dir, "src", "library_aligner", "core.py")

    cmd = [
        sys.executable, script_path,
        "-fq", str(fastq),
        "-fa", str(fasta),
        "-c", str(csv),
        "-o", str(out_bam),
        "-wc", "NNNNNNNNNNNNNNN",
        "--five_prime_flank_length", "10",
        "--three_prime_flank_length", "10",
        "--minimum_distance", "2"  # <-- CRITICAL: Now allowing up to 2 edits for the noise!
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Catch full pipeline crash
    if result.returncode != 0:
        pytest.fail(f"Pipeline Crashed!\n{result.stderr}", pytrace=False)

    pysam.index(str(out_bam))

    mapped_qnames = set()
    correct_templates = 0
    correct_strands = 0

    with pysam.AlignmentFile(out_bam, "rb") as bam:
        for read in bam.fetch():
            mapped_qnames.add(read.query_name)

            # Check if it mapped to the expected template
            if read.reference_name == "Plasmid_1":
                correct_templates += 1

            # Check if the strand orientation was handled correctly
            if "rev" in read.query_name and read.is_reverse:
                correct_strands += 1
            elif "fwd" in read.query_name and not read.is_reverse:
                correct_strands += 1

    missing_qnames = set(expected_qnames) - mapped_qnames

    # --- FAILURE REPORT ---
    if missing_qnames or correct_templates != 4 or correct_strands != 4:
        failed_reads_dump = ""
        for name, seq, _ in mappy.fastx_read(str(fastq)):
            if name in missing_qnames:
                failed_reads_dump += f">{name}\n{seq}\n\n"

        pipeline_log = [line for line in result.stderr.split('\n') if "Pipeline complete" in line]
        log_summary = pipeline_log[0] if pipeline_log else "No summary line found."

        pytest.fail(
            f"\n\n--- SLIGHT NOISE TEST FAILED ---\n"
            f"Expected 4 noisy reads to be rescued and mapped.\n"
            f"Missing reads: {missing_qnames}\n"
            f"Template Mismatches: {4 - correct_templates}\n"
            f"Strand Orientation Errors: {4 - correct_strands}\n\n"
            f"Pipeline Summary:\n{log_summary}\n\n"
            f"--- RAW SEQUENCES OF FAILED READS ---\n"
            f"Reference F5: ACGTACGTAC\n"
            f"Reference F3: TGCATGCATG\n"
            f"(Look for the 2% random mutations in the backbones/flanks below)\n\n"
            f"{failed_reads_dump}"
            f"-------------------------------------\n",
            pytrace=False
        )

    # --- SUCCESS REPORT ---
    print("\n\n" + "=" * 50)
    print("🟢 SLIGHT NOISE TEST PASSED 🟢")
    print("=" * 50)
    print(f"✅ Total noisy reads expected: {len(expected_qnames)}")
    print(f"✅ Total noisy reads rescued: {len(mapped_qnames)} / {len(expected_qnames)}")
    print(f"✅ Template accuracy: {correct_templates} / {len(mapped_qnames)} mapped to 'Plasmid_1'")
    print(f"✅ Strand awareness: {correct_strands} / {len(mapped_qnames)} oriented correctly")
    print("=" * 50 + "\n")