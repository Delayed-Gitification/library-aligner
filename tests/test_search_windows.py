import os
import subprocess
import gzip
import pysam
import pytest
import mappy


def revcomp(seq: str) -> str:
    trans = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(trans)[::-1]


def generate_window_data(tmp_path):
    """Generates 4 clean reads with barcodes heavily biased to the 5' end of the original molecule."""
    import random
    random.seed(42)  # Ensures the random DNA is the same every time you run the test

    fasta_path = tmp_path / "reference.fasta"
    csv_path = tmp_path / "linkage.csv"
    fastq_path = tmp_path / "reads.fastq.gz"
    out_bam_prefix = tmp_path / "output"

    f5 = "ACGTACGTAC"  # 10bp
    f3 = "TGCATGCATG"  # 10bp
    wildcard = "NNNNNNNNNNNNNNN"
    bases = ['A', 'C', 'G', 'T']

    # Use random complex DNA for the backbones so minimap2 anchors properly
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

        template_1 = plasmid_seq.replace(wildcard, bc_1)
        idx_1 = template_1.find(bc_1)

        # --- Shift the slice! 10bp upstream to 140bp downstream = 150bp read ---
        # Forward footprint is strictly at index [10:45]
        seq_fwd_1 = template_1[idx_1 - 10: idx_1 + 140]
        write_read("read_1_fwd_bc1", seq_fwd_1)

        seq_rev_1 = revcomp(template_1[idx_1 - 10: idx_1 + 140])
        write_read("read_2_rev_bc1", seq_rev_1)

        template_2 = plasmid_seq.replace(wildcard, bc_2)
        idx_2 = template_2.find(bc_2)

        seq_fwd_2 = template_2[idx_2 - 10: idx_2 + 140]
        write_read("read_3_fwd_bc2", seq_fwd_2)

        seq_rev_2 = revcomp(template_2[idx_2 - 10: idx_2 + 140])
        write_read("read_4_rev_bc2", seq_rev_2)

    return fasta_path, csv_path, fastq_path, out_bam_prefix, expected_qnames

def run_pipeline_with_windows(fasta, csv, fastq, out_bam, window_args):
    """Helper to run the pipeline with specific window arguments and count the results."""
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
              "--minimum_distance", "0"
          ] + window_args

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        pytest.fail(f"Pipeline Crashed!\n{result.stderr}", pytrace=False)

    pysam.index(str(out_bam))

    mapped_fwd = 0
    mapped_rev = 0
    correct_templates = 0

    with pysam.AlignmentFile(out_bam, "rb") as bam:
        for read in bam.fetch():
            if read.reference_name == "Plasmid_1":
                correct_templates += 1

            if "fwd" in read.query_name:
                mapped_fwd += 1
            elif "rev" in read.query_name:
                mapped_rev += 1

    return mapped_fwd, mapped_rev, correct_templates


def test_search_windows(tmp_path):
    fasta, csv, fastq, out_bam_prefix, expected_qnames = generate_window_data(tmp_path)

    # --- Test 1: Check FIRST 60 bases ---
    # Because of our State-Space Simplification, "first N bases" always checks the 5' end
    # of the biological molecule. It should rescue all 4 reads perfectly.
    out_1 = f"{out_bam_prefix}_first.bam"
    fwd_1, rev_1, t_1 = run_pipeline_with_windows(fasta, csv, fastq, out_1, [
        "--search_first_n_bases", "60"
    ])

    if fwd_1 != 2 or rev_1 != 2 or t_1 != 4:
        pytest.fail(f"Test 1 Failed: Expected 2 FWD and 2 REV mapped to Plasmid_1. Got {fwd_1} FWD and {rev_1} REV.",
                    pytrace=False)

    # --- Test 2: Check LAST 60 bases ---
    # The barcode is at the start of the molecule, so restricting the search to the end
    # should cause edlib to completely fail to find the flanks. All reads should drop.
    out_2 = f"{out_bam_prefix}_last.bam"
    fwd_2, rev_2, t_2 = run_pipeline_with_windows(fasta, csv, fastq, out_2, [
        "--search_last_n_bases", "60"
    ])

    if fwd_2 != 0 or rev_2 != 0:
        pytest.fail(
            f"Test 2 Failed: Expected 0 FWD and 0 REV. Got {fwd_2} FWD and {rev_2} REV. The window didn't restrict the search!",
            pytrace=False)

    # --- SUCCESS REPORT ---
    print("\n\n" + "=" * 65)
    print("🟢 NANOPORE SEARCH WINDOW BEHAVIOUR VERIFIED 🟢")
    print("=" * 65)
    print("Library Structure: Barcode footprint strictly at index [10:45]")
    print("-" * 65)
    print(f"✅ Search First 60 : Recovered {fwd_1} FWD, {rev_1} REV (Expected 2, 2)")
    print(f"✅ Search Last 60  : Recovered {fwd_2} FWD, {rev_2} REV (Expected 0, 0)")
    print("=" * 65 + "\n")