import os
import subprocess
import random
import gzip
import pysam
import pytest


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


def generate_stress_data(tmp_path, error_rate, reads_per_bc=100):
    """Generates a batch of noisy reads for a specific error rate."""
    random.seed(42)  # Keep it deterministic

    fasta_path = tmp_path / f"reference_{error_rate}.fasta"
    csv_path = tmp_path / f"linkage_{error_rate}.csv"
    fastq_path = tmp_path / f"reads_{error_rate}.fastq.gz"
    out_bam = tmp_path / f"output_{error_rate}.bam"

    bases = ['A', 'C', 'G', 'T']
    f5 = "ACGTACGTAC"
    f3 = "TGCATGCATG"
    wildcard = "NNNNNNNNNNNNNNN"

    # FIX 1: Use random complex DNA for the backbone so minimap2 anchors properly
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

    total_expected = 0

    with gzip.open(fastq_path, "wt") as f:
        def write_read(qname, seq):
            nonlocal total_expected
            total_expected += 1
            f.write(f"@{qname}\n{seq}\n+\n{'I' * len(seq)}\n")

        for bc, bc_name in [(bc_1, "bc1"), (bc_2, "bc2")]:
            template = plasmid_seq.replace(wildcard, bc)
            idx = template.find(bc)

            for i in range(reads_per_bc // 2):
                seq_fwd = template[idx - 50: idx + 100]
                noisy_fwd = add_sequencing_noise(seq_fwd, error_rate)
                write_read(f"read_{bc_name}_fwd_{i}", noisy_fwd)

                seq_rev = revcomp(template[idx - 50: idx + 100])
                noisy_rev = add_sequencing_noise(seq_rev, error_rate)
                write_read(f"read_{bc_name}_rev_{i}", noisy_rev)

    return fasta_path, csv_path, fastq_path, out_bam, total_expected


# =============================================================================
# The Automated Stress Curve Test
# =============================================================================

def test_stress_degradation_curve(tmp_path):
    import sys
    root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script_path = os.path.join(root_dir, "src", "library_aligner", "core.py")

    error_rates = [0.0, 0.02, 0.05, 0.08, 0.12, 0.16, 0.20, 0.25]
    allowed_edits = "2"

    print("\n\n" + "=" * 70)
    print(f"{'PIPELINE STRESS TEST: SENSITIVITY DEGRADATION CURVE':^70}")
    print("=" * 70)
    print(f"{'Noise %':<10} | {'Expected':<10} | {'Mapped':<10} | {'Recovery %':<12} | {'Strand OK':<10}")
    print("-" * 70)

    for rate in error_rates:
        fasta, csv, fastq, out_bam, expected_count = generate_stress_data(tmp_path, rate, reads_per_bc=100)

        cmd = [
            sys.executable, script_path,
            "-fq", str(fastq),
            "-fa", str(fasta),
            "-c", str(csv),
            "-o", str(out_bam),
            "-wc", "NNNNNNNNNNNNNNN",
            "--five_prime_flank_length", "10",
            "--three_prime_flank_length", "10",
            "--minimum_distance", allowed_edits
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            pytest.fail(f"Pipeline crashed at {rate * 100}% noise!\n{result.stderr}", pytrace=False)

        pysam.index(str(out_bam))

        # FIX 2: Use a Set to track unique reads, and filter out supplementary BAM rows
        mapped_qnames = set()
        correct_strands = 0

        with pysam.AlignmentFile(out_bam, "rb") as bam:
            for read in bam.fetch():
                # Skip secondary/supplementary hits to avoid inflating the numbers
                if read.is_secondary or read.is_supplementary:
                    continue

                mapped_qnames.add(read.query_name)

                if "rev" in read.query_name and read.is_reverse:
                    correct_strands += 1
                elif "fwd" in read.query_name and not read.is_reverse:
                    correct_strands += 1

        mapped_count = len(mapped_qnames)
        recovery_pct = (mapped_count / expected_count) * 100

        if recovery_pct >= 95:
            rec_str = f"🟢 {recovery_pct:>6.1f}%"
        elif recovery_pct >= 70:
            rec_str = f"🟡 {recovery_pct:>6.1f}%"
        else:
            rec_str = f"🔴 {recovery_pct:>6.1f}%"

        print(
            f"{rate * 100:>7.1f}%   | {expected_count:<10} | {mapped_count:<10} | {rec_str:<12} | {correct_strands}/{mapped_count}")

    print("=" * 70 + "\n")

    # We expect some loss at higher noise, but clean data should be near perfect.
    if recovery_pct < 95 and rate <= 0.02:
        pytest.fail("Sensitivity is too low on clean data!", pytrace=False)