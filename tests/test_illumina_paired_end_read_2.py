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

def generate_illumina_pe_r2_bc_data(tmp_path):
    random.seed(42)

    fasta_path = tmp_path / "reference_pe_r2bc.fasta"
    csv_path = tmp_path / "linkage_pe_r2bc.csv"
    fastq_r1_path = tmp_path / "reads_R1.fastq.gz"
    fastq_r2_path = tmp_path / "reads_R2.fastq.gz"
    out_bam = tmp_path / "output_pe_r2bc.bam"

    bases = ['A', 'C', 'G', 'T']
    ref_A = "".join(random.choices(bases, k=600))
    ref_B = "".join(random.choices(bases, k=600))

    with open(fasta_path, "w") as f:
        f.write(f">Genome_A\n{ref_A}\n")
        f.write(f">Genome_B\n{ref_B}\n")

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
    # R1: 150bp pure genomic (forward strand, from the start of the insert)
    # R2: 10bp adapter + 15bp barcode + 125bp genomic (reverse complemented insert end)
    #
    # This mirrors a library prep where the barcode is ligated to the R2 end,
    # e.g. a 3' tag-seq or a UMI/barcode on the second adapter.
    adapter = "AATGATACGGCGACC"  # 15bp — barcode lands at R2[15:30]
    insert_size = 250
    read_len = 150
    genomic_in_r2 = read_len - len(adapter) - len(bc_A1)  # = 125bp

    bc_start = len(adapter)             # 15
    bc_end = len(adapter) + len(bc_A1)  # 30

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

                max_start = len(ref_seq) - insert_size
                start_idx = random.randint(0, max_start)
                insert = ref_seq[start_idx: start_idx + insert_size]

                # R1: first 150bp of insert, forward strand, no adapter/barcode
                r1_genomic = insert[:read_len]
                r1 = r1_genomic

                # R2: adapter + barcode + last 125bp of insert, reverse complemented
                # The sequencer reads the opposite strand from the other end,
                # so the genomic portion is revcomp'd, but the adapter+barcode
                # are prepended in the forward orientation as the sequencer sees them
                r2_genomic = revcomp(insert[insert_size - genomic_in_r2:])
                r2 = adapter + bc + r2_genomic

                r1_noisy = add_illumina_noise(r1, error_rate=0.005)
                r2_noisy = add_illumina_noise(r2, error_rate=0.005)

                write_pair(f"read_{read_counter}_EXPECT_{rname}", r1_noisy, r2_noisy, rname)

    return fasta_path, csv_path, fastq_r1_path, fastq_r2_path, out_bam, expected_qnames, ground_truth, bc_start, bc_end


# =============================================================================
# The Test
# =============================================================================

def test_illumina_paired_end_r2_barcode(tmp_path):
    fasta, csv, fastq_r1, fastq_r2, out_bam, expected_qnames, ground_truth, bc_start, bc_end = generate_illumina_pe_r2_bc_data(tmp_path)

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
        "--barcode_read", "2",               # <-- key difference
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

    complete_pairs = r1_seen & r2_seen
    missing_qnames = set(expected_qnames) - mapped_qnames
    mapped_val = len(mapped_qnames)
    target_val = len(expected_qnames) * 0.95

    if mapped_val < target_val or correct_templates < mapped_val * 2 * 0.95:
        failed_reads_dump = ""
        dump_count = 0
        for name, seq, _ in mappy.fastx_read(str(fastq_r2)):
            if name in missing_qnames and dump_count < 5:
                failed_reads_dump += f">{name} (R2)\n{seq}\n\n"
                dump_count += 1

        pipeline_log = [line for line in result.stderr.split('\n') if "Pipeline complete" in line]
        log_summary = pipeline_log[0] if pipeline_log else "No summary line found."

        pytest.fail(
            f"\n\n--- ILLUMINA PAIRED-END (R2 BARCODE) TEST FAILED ---\n"
            f"Expected {len(expected_qnames)} read pairs. Recovered: {mapped_val}/{len(expected_qnames)}\n"
            f"Complete pairs (R1+R2 both present): {len(complete_pairs)}/{len(expected_qnames)}\n"
            f"Template mismatches: {mapped_val * 2 - correct_templates}\n\n"
            f"Pipeline Summary:\n{log_summary}\n\n"
            f"--- SAMPLE OF FAILED R2 READS ---\n"
            f"{failed_reads_dump}"
            f"----------------------------------------------------\n",
            pytrace=False
        )

    print("\n\n" + "=" * 70)
    print("🟢 ILLUMINA PAIRED-END (R2 BARCODE) TEST PASSED 🟢")
    print("=" * 70)
    print(f"✅ Configuration: PE150 | Barcode on R2[{bc_start}:{bc_end}] | Preset: sr")
    print(f"✅ Total pairs expected:       {len(expected_qnames)}")
    print(f"✅ Total pairs recovered:      {mapped_val} / {len(expected_qnames)}")
    print(f"✅ Complete pairs (R1+R2):     {len(complete_pairs)} / {len(expected_qnames)}")
    print(f"✅ Template accuracy:          {correct_templates} / {mapped_val * 2} alignments")
    print("=" * 70 + "\n")