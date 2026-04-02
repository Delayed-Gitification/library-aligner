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


def generate_absent_barcode_data(tmp_path):
    """
    Generates Option 2 (CSV) data:
    Barcodes map to templates via CSV, but the FASTA templates DO NOT contain the barcode.
    """
    random.seed(42)

    fasta_path = tmp_path / "reference.fasta"
    csv_path = tmp_path / "linkage.csv"
    fastq_path = tmp_path / "reads.fastq.gz"
    out_bam = tmp_path / "output.bam"

    bases = ['A', 'C', 'G', 'T']
    f5 = "ACGTACGTAC"  # 10bp
    f3 = "TGCATGCATG"  # 10bp

    # 1. Create the References (NO BARCODES, NO FLANKS, JUST BACKBONE)
    # We use random DNA to ensure minimap2 anchors properly
    backbone_A_L = "".join(random.choices(bases, k=150))
    backbone_A_R = "".join(random.choices(bases, k=150))
    ref_A = backbone_A_L + backbone_A_R

    backbone_B_L = "".join(random.choices(bases, k=150))
    backbone_B_R = "".join(random.choices(bases, k=150))
    ref_B = backbone_B_L + backbone_B_R

    with open(fasta_path, "w") as f:
        f.write(f">Plasmid_A\n{ref_A}\n")
        f.write(f">Plasmid_B\n{ref_B}\n")

    # 2. Assign multiple barcodes per template
    bc_A1 = "GATTACAGATTACAG"
    bc_A2 = "CCATCCATCCATCCA"
    bc_B1 = "TTAGCTTAGCTTAGC"
    bc_B2 = "GGAAGGAAGGAAGGA"

    with open(csv_path, "w") as f:
        f.write("barcode,rname\n")
        f.write(f"{bc_A1},Plasmid_A\n")
        f.write(f"{bc_A2},Plasmid_A\n")
        f.write(f"{bc_B1},Plasmid_B\n")
        f.write(f"{bc_B2},Plasmid_B\n")

    expected_qnames = []

    # We will track exactly what each read is supposed to map to
    ground_truth = {}

    with gzip.open(fastq_path, "wt") as f:
        def write_read(qname, seq, expected_rname, strand):
            expected_qnames.append(qname)
            ground_truth[qname] = {"rname": expected_rname, "strand": strand}
            f.write(f"@{qname}\n{seq}\n+\n{'I' * len(seq)}\n")

        # 3. Generate Reads: Backbone + Flank + Barcode + Flank + Backbone
        # Even though the reference lacks the flanks/bc, the reads MUST have them
        configs = [
            (bc_A1, "Plasmid_A", backbone_A_L, backbone_A_R),
            (bc_A2, "Plasmid_A", backbone_A_L, backbone_A_R),
            (bc_B1, "Plasmid_B", backbone_B_L, backbone_B_R),
            (bc_B2, "Plasmid_B", backbone_B_L, backbone_B_R)
        ]

        read_counter = 0
        reads_per_bc = 100  # 50 fwd, 50 rev per barcode = 400 total reads

        for bc, rname, left_bb, right_bb in configs:
            # The true biological molecule
            true_molecule = left_bb + f5 + bc + f3 + right_bb

            for _ in range(reads_per_bc // 2):
                read_counter += 1

                # Forward Read
                noisy_fwd = add_sequencing_noise(true_molecule, error_rate=0.02)
                write_read(f"read_{read_counter}_EXPECT_{rname}_STRAND_fwd", noisy_fwd, rname, "fwd")

                # Reverse Read
                read_counter += 1
                noisy_rev = add_sequencing_noise(revcomp(true_molecule), error_rate=0.02)
                write_read(f"read_{read_counter}_EXPECT_{rname}_STRAND_rev", noisy_rev, rname, "rev")

    return fasta_path, csv_path, fastq_path, out_bam, expected_qnames, ground_truth


# =============================================================================
# The Diagnostic Test
# =============================================================================

def test_option2_csv_absent_barcodes(tmp_path):
    fasta, csv, fastq, out_bam, expected_qnames, ground_truth = generate_absent_barcode_data(tmp_path)

    import sys
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script_path = os.path.join(root_dir, "src", "library_aligner", "core.py")

    cmd = [
        sys.executable, script_path,
        "-fq", str(fastq),
        "-fa", str(fasta),
        "-c", str(csv),  # <-- CRITICAL: Using CSV for linkage
        "-o", str(out_bam),
        # NO WILDCARD TAG (-wc) INCLUDED HERE!
        "--five_prime_flank", "ACGTACGTAC",  # <-- CRITICAL: Explicit flanks to find BCs
        "--three_prime_flank", "TGCATGCATG",
        "--minimum_distance", "2"
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
            # Filter out secondary/supplementary alignments
            if read.is_secondary or read.is_supplementary:
                continue

            qname = read.query_name
            mapped_qnames.add(qname)

            # Fetch the ground truth for this specific read
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
    target_val = len(expected_qnames) * 0.95  # Require 95% recovery

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
            f"\n\n--- OPTION 2 (ABSENT BARCODES) TEST FAILED ---\n"
            f"Expected 400 noisy reads. Recovered: {mapped_val}/400\n"
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
    print("🟢 OPTION 2 (CSV + ABSENT BARCODES) TEST PASSED 🟢")
    print("=" * 70)
    print(f"✅ Configuration: FASTA has no barcodes. CSV links multiple BCs per FASTA.")
    print(f"✅ Total noisy reads expected: {len(expected_qnames)}")
    print(f"✅ Total noisy reads rescued: {mapped_val} / {len(expected_qnames)} (>95%)")
    print(f"✅ Multi-Mapping accuracy: {correct_templates} / {mapped_val} mapped to correct template")
    print(f"✅ Strand awareness: {correct_strands} / {mapped_val} oriented correctly")
    print("=" * 70 + "\n")