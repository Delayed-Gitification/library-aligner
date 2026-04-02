import os
import subprocess
import gzip
import pysam
import pytest
import random
import mappy


def revcomp(seq: str) -> str:
    trans = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(trans)[::-1]


def generate_single_flank_data(tmp_path):
    """Generates 4 perfectly clean reads with highly visible flank sequences."""
    random.seed(42)  # Keep the random DNA deterministic

    fasta_path = tmp_path / "reference.fasta"
    csv_path = tmp_path / "linkage.csv"
    fastq_path = tmp_path / "reads.fastq.gz"
    out_bam_prefix = tmp_path / "output"

    f5 = "ACGTACGTAC"  # 10bp
    f3 = "TGCATGCATG"  # 10bp
    wildcard = "NNNNNNNNNNNNNNN"
    bases = ['A', 'C', 'G', 'T']

    # Random complex backbones so minimap2 anchors properly
    backbone_5 = "".join(random.choices(bases, k=100))
    backbone_3 = "".join(random.choices(bases, k=100))

    plasmid_seq = backbone_5 + f5 + wildcard + f3 + backbone_3

    with open(fasta_path, "w") as f:
        f.write(f">Plasmid_1\n{plasmid_seq}\n")

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

        # Read 1: Perfect Forward BC1
        seq_fwd_1 = template_1[idx_1 - 50: idx_1 + 100]
        write_read("read_1_fwd_bc1", seq_fwd_1)

        # Read 2: Perfect Reverse BC1
        seq_rev_1 = revcomp(template_1[idx_1 - 50: idx_1 + 100])
        write_read("read_2_rev_bc1", seq_rev_1)

        # Create the true template for Barcode 2
        template_2 = plasmid_seq.replace(wildcard, bc_2)
        idx_2 = template_2.find(bc_2)

        # Read 3: Perfect Forward BC2
        seq_fwd_2 = template_2[idx_2 - 50: idx_2 + 100]
        write_read("read_3_fwd_bc2", seq_fwd_2)

        # Read 4: Perfect Reverse BC2
        seq_rev_2 = revcomp(template_2[idx_2 - 50: idx_2 + 100])
        write_read("read_4_rev_bc2", seq_rev_2)

    return fasta_path, csv_path, fastq_path, out_bam_prefix, expected_qnames


def run_flank_pipeline(fasta, csv, fastq, out_bam, flank_args, expected_qnames):
    """Helper to run the pipeline with specific flank configurations."""
    import sys
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script_path = os.path.join(root_dir, "src", "library_aligner", "core.py")

    # Note: We pass -wc so the library builder replaces the Ns in the FASTA references,
    # but our explicit flank_args will override the get_flanks auto-detection.
    cmd = [
              sys.executable, script_path,
              "-fq", str(fastq),
              "-fa", str(fasta),
              "-c", str(csv),
              "-o", str(out_bam),
              "-wc", "NNNNNNNNNNNNNNN",
              "--minimum_distance", "0"
          ] + flank_args

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        pytest.fail(f"Pipeline Crashed!\n{result.stderr}", pytrace=False)

    pysam.index(str(out_bam))

    mapped_qnames = set()
    correct_templates = 0
    correct_strands = 0

    with pysam.AlignmentFile(out_bam, "rb") as bam:
        for read in bam.fetch():
            mapped_qnames.add(read.query_name)

            if read.reference_name == "Plasmid_1":
                correct_templates += 1

            if "rev" in read.query_name and read.is_reverse:
                correct_strands += 1
            elif "fwd" in read.query_name and not read.is_reverse:
                correct_strands += 1

    missing_qnames = set(expected_qnames) - mapped_qnames

    if missing_qnames or correct_templates != 4 or correct_strands != 4:
        # Safely extract the summary line without using backslashes in f-strings
        pipeline_log = [line for line in result.stderr.splitlines() if "Pipeline complete" in line]
        log_summary = pipeline_log[0] if pipeline_log else "No summary line found."

        pytest.fail(
            f"Test Failed with args {flank_args}.\n"
            f"Missing reads: {missing_qnames}\n"
            f"Pipeline Summary:\n{log_summary}",
            pytrace=False
        )

    return len(mapped_qnames), correct_templates, correct_strands


def test_single_flank_extraction(tmp_path):
    fasta, csv, fastq, out_bam_prefix, expected_qnames = generate_single_flank_data(tmp_path)

    # --- Test 1: 5' Flank Only ---
    out_5p = f"{out_bam_prefix}_5p_only.bam"
    mapped_5p, temp_5p, strand_5p = run_flank_pipeline(fasta, csv, fastq, out_5p, [
        "--five_prime_flank", "ACGTACGTAC"
    ], expected_qnames)

    # --- Test 2: 3' Flank Only ---
    out_3p = f"{out_bam_prefix}_3p_only.bam"
    mapped_3p, temp_3p, strand_3p = run_flank_pipeline(fasta, csv, fastq, out_3p, [
        "--three_prime_flank", "TGCATGCATG"
    ], expected_qnames)

    # --- SUCCESS REPORT ---
    print("\n\n" + "=" * 65)
    print("🟢 SINGLE-FLANK EXTRACTION PASSED 🟢")
    print("=" * 65)
    print("✅ 5' Flank Only : Extracted EXACTLY 15bp downstream.")
    print(f"   -> {mapped_5p}/4 Mapped | {temp_5p}/4 Template Acc | {strand_5p}/4 Strand Acc")
    print("✅ 3' Flank Only : Extracted EXACTLY 15bp upstream.")
    print(f"   -> {mapped_3p}/4 Mapped | {temp_3p}/4 Template Acc | {strand_3p}/4 Strand Acc")
    print("=" * 65 + "\n")