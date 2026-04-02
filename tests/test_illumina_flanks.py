import os
import subprocess
import random
import gzip
import pysam
import pytest
import mappy


# =============================================================================
# Helpers
# =============================================================================

def revcomp(seq: str) -> str:
    trans = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(trans)[::-1]


def add_illumina_noise(seq: str, error_rate: float) -> str:
    bases = ['A', 'C', 'G', 'T']
    noisy_seq = []
    for base in seq:
        if random.random() < error_rate:
            noisy_seq.append(random.choice([b for b in bases if b != base]))
        else:
            noisy_seq.append(base)
    return "".join(noisy_seq)


# =============================================================================
# Data Generation
# =============================================================================

def generate_pe_flank_data(tmp_path):
    random.seed(42)

    fasta_path = tmp_path / "reference_pe_flank.fasta"
    csv_path = tmp_path / "linkage_pe_flank.csv"
    fastq_r1_path = tmp_path / "reads_R1.fastq.gz"
    fastq_r2_path = tmp_path / "reads_R2.fastq.gz"
    out_bam = tmp_path / "output_pe_flank.bam"

    # Flanks and wildcard — identical structure to your single-end flank test
    f5 = "ACGTACGTAC"   # 10bp
    f3 = "TGCATGCATG"   # 10bp
    wildcard = "NNNNNNNNNNNNNNN"  # 15bp placeholder

    # Plasmid backbone: fixed regions either side of the barcode cassette.
    # Made long enough that PE150 reads have plenty of genomic content to align.
    bases = ['A', 'C', 'G', 'T']
    left_backbone  = "".join(random.choices(bases, k=200))
    right_backbone = "".join(random.choices(bases, k=200))

    # Template with wildcard — written to FASTA as-is for the pipeline
    plasmid_template = left_backbone + f5 + wildcard + f3 + right_backbone

    bc_1 = "GATTACAGATTACAG"
    bc_2 = "CCATCCATCCATCCA"
    bc_3 = "TTAGCTTAGCTTAGC"
    bc_4 = "GGAAGGAAGGAAGGA"

    with open(fasta_path, "w") as f:
        f.write(f">Plasmid_1\n{plasmid_template}\n")

    with open(csv_path, "w") as f:
        f.write("barcode,rname\n")
        for bc in [bc_1, bc_2, bc_3, bc_4]:
            f.write(f"{bc},Plasmid_1\n")

    expected_qnames = []
    ground_truth = {}

    # For each barcode, build the resolved plasmid sequence (wildcard replaced)
    # R1: 150bp window centred on the barcode cassette (adapter + f5 + bc + f3 + some backbone)
    # R2: 150bp from the right backbone, reverse complemented — pure genomic, no barcode
    read_len = 150
    insert_size = 300

    with gzip.open(fastq_r1_path, "wt") as f1, gzip.open(fastq_r2_path, "wt") as f2:

        def write_pair(qname, r1_seq, r2_seq):
            expected_qnames.append(qname)
            ground_truth[qname] = {"rname": "Plasmid_1"}
            f1.write(f"@{qname}\n{r1_seq}\n+\n{'I' * len(r1_seq)}\n")
            f2.write(f"@{qname}\n{r2_seq}\n+\n{'I' * len(r2_seq)}\n")

        read_counter = 0
        reads_per_bc = 100

        for bc in [bc_1, bc_2, bc_3, bc_4]:
            # Build the fully resolved plasmid for this barcode
            resolved = plasmid_template.replace(wildcard, bc)

            # Find the barcode position in the resolved sequence
            bc_idx = resolved.find(bc)

            # R1 starts 50bp upstream of the barcode so the flanks are well-centred
            # within the read — mirrors your single-end flank test geometry
            r1_start = bc_idx - 50
            r1_end = r1_start + read_len

            # R2 comes from insert_size bp downstream of r1_start, revcomp'd
            r2_end = r1_start + insert_size
            r2_start = r2_end - read_len

            for i in range(reads_per_bc):
                read_counter += 1

                r1_seq = resolved[r1_start:r1_end]
                r2_seq = revcomp(resolved[r2_start:r2_end])

                r1_noisy = add_illumina_noise(r1_seq, error_rate=0.01)
                r2_noisy = add_illumina_noise(r2_seq, error_rate=0.01)

                write_pair(f"read_{read_counter}_EXPECT_Plasmid_1", r1_noisy, r2_noisy)

    return fasta_path, csv_path, fastq_r1_path, fastq_r2_path, out_bam, expected_qnames, ground_truth


# =============================================================================
# The Test
# =============================================================================

def test_pe_flank_based_barcode(tmp_path):
    fasta, csv, fastq_r1, fastq_r2, out_bam, expected_qnames, ground_truth = generate_pe_flank_data(tmp_path)

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
        "-wc", "NNNNNNNNNNNNNNN",
        "--five_prime_flank_length", "10",
        "--three_prime_flank_length", "10",
        "--minimum_distance", "1",
        "--minimap2_preset", "sr",
        # No --barcode_start_pos, no --barcode_read — flank-based path, R1 by default
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

            if read.reference_name == ground_truth[qname]["rname"]:
                correct_templates += 1

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
            f"\n\n--- PE FLANK-BASED BARCODE TEST FAILED ---\n"
            f"Expected {len(expected_qnames)} read pairs. Recovered: {mapped_val}/{len(expected_qnames)}\n"
            f"Complete pairs (R1+R2 both present): {len(complete_pairs)}/{len(expected_qnames)}\n"
            f"Template mismatches: {mapped_val * 2 - correct_templates}\n\n"
            f"Pipeline Summary:\n{log_summary}\n\n"
            f"--- SAMPLE OF FAILED R1 READS ---\n"
            f"{failed_reads_dump}"
            f"------------------------------------------\n",
            pytrace=False
        )

    print("\n\n" + "=" * 70)
    print("🟢 PE FLANK-BASED BARCODE TEST PASSED 🟢")
    print("=" * 70)
    print(f"✅ Configuration: PE150 | Flank-based extraction | Wildcard: NNNNNNNNNNNNNNN")
    print(f"✅ Flanks: 5'=ACGTACGTAC  3'=TGCATGCATG | Max edits: 2")
    print(f"✅ Total pairs expected:       {len(expected_qnames)}")
    print(f"✅ Total pairs recovered:      {mapped_val} / {len(expected_qnames)}")
    print(f"✅ Complete pairs (R1+R2):     {len(complete_pairs)} / {len(expected_qnames)}")
    print(f"✅ Template accuracy:          {correct_templates} / {mapped_val * 2} alignments")
    print("=" * 70 + "\n")