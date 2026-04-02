import os
import subprocess
import gzip
import pysam
import pytest
import mappy
import random


def revcomp(seq: str) -> str:
    trans = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(trans)[::-1]


# =============================================================================
# Spliced Data Generation
# =============================================================================

def generate_spliced_data(tmp_path):
    """
    Generates a reference with a 100nt intron (GT-AG) and
    RNA reads where that intron has been spliced out.
    """
    random.seed(42)
    fasta_path = tmp_path / "reference_spliced.fasta"
    csv_path = tmp_path / "linkage_spliced.csv"
    fastq_path = tmp_path / "reads_spliced.fastq.gz"
    out_bam = tmp_path / "output_spliced.bam"

    f5 = "ACGTACGTAC"  # 10bp
    f3 = "TGCATGCATG"  # 10bp
    wildcard = "NNNNNNNNNNNNNNN"  # 15bp

    # Build Exons and Intron
    # Exon 1 contains the barcode cassette
    exon1_seq = ("A" * 50) + f5 + wildcard + f3 + ("A" * 50)  # 135bp

    # 100nt Intron with canonical splice sites
    intron_seq = "GT" + "".join(random.choices("ACGT", k=96)) + "AG"  # 100bp

    exon2_seq = ("T" * 100)  # 100bp

    # Full genomic reference sequence
    full_ref = exon1_seq + intron_seq + exon2_seq

    with open(fasta_path, "w") as f:
        f.write(f">Spliced_Ref\n{full_ref}\n")

    bc_1 = "GATTACAGATTACAG"
    with open(csv_path, "w") as f:
        f.write("barcode,rname\n")
        f.write(f"{bc_1},Spliced_Ref\n")

    expected_qnames = []
    with gzip.open(fastq_path, "wt") as f:
        def write_read(qname, seq):
            expected_qnames.append(qname)
            f.write(f"@{qname}\n{seq}\n+\n{'I' * len(seq)}\n")

        # The "Mature mRNA" read (Exon 1 + Exon 2)
        mature_template = (exon1_seq + exon2_seq).replace(wildcard, bc_1)

        # Read 1: Forward mature mRNA (full length)
        write_read("read_fwd_spliced", mature_template)

        # Read 2: Reverse mature mRNA
        write_read("read_rev_spliced", revcomp(mature_template))

    return fasta_path, csv_path, fastq_path, out_bam, expected_qnames


# =============================================================================
# Test Logic
# =============================================================================

def test_spliced_alignment_capability(tmp_path):
    fasta, csv, fastq, out_bam, expected_qnames = generate_spliced_data(tmp_path)

    import sys
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script_path = os.path.join(root_dir, "src", "library_aligner", "core.py")

    # CRITICAL: We MUST use -ax splice or ensure the preset allows for gaps
    cmd = [
        sys.executable, script_path,
        "-fq", str(fastq),
        "-fa", str(fasta),
        "-c", str(csv),
        "-o", str(out_bam),
        "-wc", "NNNNNNNNNNNNNNN",
        "--five_prime_flank_length", "10",
        "--three_prime_flank_length", "10",
        "--minimap2_preset", "splice",  # Ensure splice mode is active
        "--minimum_distance", "0"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        pytest.fail(f"Pipeline Crashed!\n{result.stderr}", pytrace=False)

    pysam.index(str(out_bam))

    spliced_reads_found = 0
    correct_intron_size = 0

    with pysam.AlignmentFile(out_bam, "rb") as bam:
        for read in bam.fetch():
            print("")
            print(read)
            # CIGAR op 3 is 'N' (ref_skip), which represents the intron
            cigar = read.cigartuples
            if cigar:
                for op, length in cigar:
                    if op == 3:  # Found a splice junction
                        spliced_reads_found += 1
                        if length == 100:
                            correct_intron_size += 1

    # --- FAILURE ANALYSIS ---
    if spliced_reads_found < 2 or correct_intron_size < 2:
        pytest.fail(
            f"\n\n--- SPLICED ALIGNMENT TEST FAILED ---\n"
            f"Expected 2 spliced alignments with 100nt gaps.\n"
            f"Spliced alignments found: {spliced_reads_found}\n"
            f"Correct intron sizes (100nt): {correct_intron_size}\n"
            f"Check if --minimap2_preset is correctly passing 'splice'.\n"
            f"----------------------------------------\n",
            pytrace=False
        )

    print("\n\n" + "=" * 50)
    print("🟢 SPLICED ALIGNMENT TEST PASSED 🟢")
    print("=" * 50)
    print(f"✅ Introns detected via CIGAR 'N' ops: {spliced_reads_found}")
    print(f"✅ Accurate 100nt intron recovery:    {correct_intron_size}")
    print("=" * 50 + "\n")