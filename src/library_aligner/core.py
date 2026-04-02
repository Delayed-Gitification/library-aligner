import os
import tempfile
import logging
import bisect
from dataclasses import dataclass
from enum import Enum
from typing import Optional
import mappy
import edlib
import pysam
import pandas as pd
import argparse

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


# =============================================================================
# Data Models & Enums
# =============================================================================

class MappingResult(Enum):
    SUCCESS = "success"
    NO_FLANKS_FOUND = "no_flanks_found"
    BARCODE_UNRECOGNISED = "barcode_unrecognised"
    BARCODE_AMBIGUOUS = "barcode_ambiguous"
    MAPPING_FAILED = "mapping_failed"


@dataclass
class PipelineStats:
    processed: int = 0
    no_flanks_found: int = 0
    barcode_unrecognised: int = 0
    barcode_ambiguous: int = 0
    mapping_failed: int = 0
    mapped: int = 0

    def log_summary(self):
        logger.info(
            f"Pipeline complete: {self.processed:,} reads processed | "
            f"{self.mapped:,} mapped | "
            f"{self.no_flanks_found:,} no flanks | "
            f"{self.barcode_unrecognised:,} barcode unrecognised | "
            f"{self.barcode_ambiguous:,} barcode ambiguous | "
            f"{self.mapping_failed:,} mapping failed"
        )


# =============================================================================
# Core Classes & Logic
# =============================================================================

class BarcodeMatcher:
    def __init__(self, known_barcodes: list[str]):
        self.known_barcodes = known_barcodes
        self.exact_set = set(known_barcodes)

    def match(self, extracted_seq: str, max_edits: int = 2) -> tuple[Optional[str], MappingResult]:
        if extracted_seq in self.exact_set:
            return extracted_seq, MappingResult.SUCCESS

        if max_edits < 1:
            return None, MappingResult.BARCODE_UNRECOGNISED

        best_matches = set()
        # Fallback to standard 1-to-1 comparison for fuzzy reads
        for bc in self.known_barcodes:
            # "NW" means global alignment: it must match end-to-end
            res = edlib.align(extracted_seq, bc, mode="NW", k=max_edits)
            if res["editDistance"] != -1:
                best_matches.add(bc)

        if not best_matches:
            return None, MappingResult.BARCODE_UNRECOGNISED

        if len(best_matches) > 1:
            return None, MappingResult.BARCODE_AMBIGUOUS

        return best_matches.pop(), MappingResult.SUCCESS


def extract_barcode_sequence(
        read_seq: str, f5_fwd=None, f3_fwd=None, f5_rev=None, f3_rev=None,
        expected_bc_len=15, search_window=350, max_error_rate=0.2,
        start_pos=None, end_pos=None,
) -> Optional[str]:
    if start_pos is not None and end_pos is not None:
        return read_seq[start_pos:end_pos] if end_pos <= len(read_seq) else None

    if not (f5_fwd and f3_fwd): return None

    seq_len = len(read_seq)
    search_window = min(search_window, seq_len)
    windows = [read_seq[:search_window], read_seq[-search_window:]] if seq_len > search_window else [read_seq]

    fwd_template = f5_fwd + ("N" * expected_bc_len) + f3_fwd
    rev_template = f5_rev + ("N" * expected_bc_len) + f3_rev
    max_errs = int(len(fwd_template) * max_error_rate)
    dna_eq = [("N", "A"), ("N", "C"), ("N", "G"), ("N", "T")]

    for window_seq in windows:
        for orient, template, f5_s, f3_s in [("f", fwd_template, f5_fwd, f3_fwd), ("r", rev_template, f5_rev, f3_rev)]:
            res = edlib.align(template, window_seq, "HW", "locations", max_errs, additionalEqualities=dna_eq)
            if res["editDistance"] == -1: continue

            s, e = res["locations"][0]
            footprint = window_seq[max(0, s - 3): min(len(window_seq), e + 4)]

            # FIX: We must give edlib tolerance for the flanks so noisy reads survive!
            f5_errs = int(len(f5_s) * max_error_rate)
            f3_errs = int(len(f3_s) * max_error_rate)

            res_5 = edlib.align(f5_s, footprint, "HW", "locations", k=f5_errs)
            res_3 = edlib.align(f3_s, footprint, "HW", "locations", k=f3_errs)

            if res_5["editDistance"] != -1 and res_3["editDistance"] != -1:
                e5, s3 = res_5["locations"][0][1], res_3["locations"][-1][0]
                if s3 > e5:
                    extracted = footprint[e5 + 1: s3]
                    return mappy.revcomp(extracted) if orient == "r" else extracted
    return None

# =============================================================================
# BAM Handling
# =============================================================================

def _build_aligned_segment(
        bam_writer: pysam.AlignmentFile,
        ref_id: int,
        read_name: str,
        read_seq: str,
        aln: mappy.Alignment,
        cached_qual_array,
) -> pysam.AlignedSegment:
    """Constructs a pysam AlignedSegment using SOFT clipping."""
    a = pysam.AlignedSegment(bam_writer.header)
    a.query_name = read_name
    a.reference_id = ref_id
    a.reference_start = aln.r_st
    a.mapping_quality = aln.mapq
    a.is_supplementary = not aln.is_primary
    a.is_reverse = (aln.strand == -1)

    read_len = len(read_seq)
    clip_left = aln.q_st
    clip_right = read_len - aln.q_en

    cigar = [(op, length) for length, op in aln.cigar]

    if aln.strand == 1:
        if clip_left > 0:
            cigar.insert(0, (4, clip_left))
        if clip_right > 0:
            cigar.append((4, clip_right))
    else:
        if clip_right > 0:
            cigar.insert(0, (4, clip_right))
        if clip_left > 0:
            cigar.append((4, clip_left))

    a.cigar = cigar

    if a.is_supplementary:
        a.query_sequence = None
        a.query_qualities = None
    else:
        if a.is_reverse:
            a.query_sequence = mappy.revcomp(read_seq)
            if cached_qual_array is not None:
                a.query_qualities = cached_qual_array[::-1]
        else:
            a.query_sequence = read_seq
            if cached_qual_array is not None:
                a.query_qualities = cached_qual_array
    return a


def process_and_write_read(bam_writer, ref_id, name, seq, qual, barcode, aligners_dict,
                           min_mapq=0):  # FIX: Dropped min_mapq to 0 for mini-reads
    aligner = aligners_dict.get(barcode)
    if not aligner: return MappingResult.BARCODE_UNRECOGNISED

    alns = [a for a in aligner.map(seq) if a.mapq >= min_mapq]
    if not alns: return MappingResult.MAPPING_FAILED

    q_arr = pysam.qualitystring_to_array(qual) if qual else None
    for aln in alns:
        a = pysam.AlignedSegment(bam_writer.header)
        a.query_name, a.reference_id, a.reference_start, a.mapping_quality = name, ref_id, aln.r_st, aln.mapq
        a.is_reverse = (aln.strand == -1)

        cigar = [(op, length) for length, op in aln.cigar]
        clip_l, clip_r = aln.q_st, len(seq) - aln.q_en
        if aln.strand == 1:
            if clip_l > 0: cigar.insert(0, (4, clip_l))
            if clip_r > 0: cigar.append((4, clip_r))
        else:
            if clip_r > 0: cigar.insert(0, (4, clip_r))
            if clip_l > 0: cigar.append((4, clip_l))
        a.cigar = cigar

        if not (not aln.is_primary):
            if a.is_reverse:
                a.query_sequence = mappy.revcomp(seq)
                if q_arr is not None: a.query_qualities = q_arr[::-1]
            else:
                a.query_sequence = seq
                if q_arr is not None: a.query_qualities = q_arr
        bam_writer.write(a)
    return MappingResult.SUCCESS


# =============================================================================
# Setup Helpers
# =============================================================================

def build_aligner_dict(
        plasmid_references: dict[str, dict[str, str]],
        preset: str = "splice",
        kmer: Optional[int] = None,
        w_score: Optional[int] = None
) -> dict[str, mappy.Aligner]:
    """Initialises splice-aware mappy aligners using custom settings."""
    aligners = {}

    # Clean up common CLI prefixes if the user types "-ax splice" instead of "splice"
    clean_preset = preset.replace("-ax ", "").replace("-x ", "").strip()

    for barcode, ref_data in plasmid_references.items():
        # Build the base arguments
        kwargs = {
            "seq": ref_data["seq"],
            "preset": clean_preset,
            "min_chain_score": 25
        }

        # Only add k and w if the user explicitly provided them
        if kmer is not None:
            kwargs["k"] = kmer
        if w_score is not None:
            kwargs["w"] = w_score

        # Unpack the dictionary into the Aligner
        aligner = mappy.Aligner(**kwargs)

        if not aligner:
            raise RuntimeError(f"Failed to load reference for barcode: {barcode} using preset '{clean_preset}'")

        aligners[barcode] = aligner

    return aligners


def create_bam_header(plasmid_references: dict[str, dict[str, str]]) -> tuple[dict, dict[str, int]]:
    """Builds unified SAM header with unique reference names."""
    header = {"HD": {"VN": "1.0", "SO": "unsorted"}, "SQ": []}
    ref_name_to_id: dict[str, int] = {}
    unique_refs = {}

    for ref_data in plasmid_references.values():
        rname = ref_data["rname"]
        if rname not in unique_refs:
            unique_refs[rname] = len(ref_data["seq"])

    for i, (rname, ref_len) in enumerate(sorted(unique_refs.items())):
        header["SQ"].append({"SN": rname, "LN": ref_len})
        ref_name_to_id[rname] = i

    return header, ref_name_to_id


def build_library_dictionary(
        csv_path: Optional[str],
        fasta_path: str,
        to_replace: Optional[str] = None,
        rname_equals_barcode: bool = False
) -> dict[str, dict[str, str]]:
    """Links barcodes to their reference sequences, either via CSV or direct FASTA headers."""

    plasmid_library = {}

    # Option 2 Fast-Path: FASTA headers *are* the barcodes
    if rname_equals_barcode:
        logger.info("Using FASTA headers directly as barcodes (--rname_equals_barcode).")
        for name, base_seq, _ in mappy.fastx_read(fasta_path):
            # The FASTA header is the barcode, and there is no wildcard replacement needed
            plasmid_library[name] = {"rname": name, "seq": base_seq}
        return plasmid_library

    # Option 1: Standard CSV Linkage
    if not csv_path:
        raise ValueError("A CSV file must be provided unless --rname_equals_barcode is active.")

    logger.info("Parsing CSV to link barcodes to FASTA references...")
    df = pd.read_csv(csv_path)

    rname_to_barcodes = {}
    for _, row in df.iterrows():
        rname_to_barcodes.setdefault(row.rname, []).append(row.barcode)

    for name, base_seq, _ in mappy.fastx_read(fasta_path):
        if name in rname_to_barcodes:
            for barcode in rname_to_barcodes[name]:
                # Inject the barcode into the wildcard position
                seq = base_seq.replace(to_replace, barcode) if to_replace else base_seq
                plasmid_library[barcode] = {"rname": name, "seq": seq}
        else:
            logger.warning(f"Reference '{name}' found in FASTA but missing from CSV linkage.")

    return plasmid_library


def get_flanks(
        fasta_path: str,
        bc_placeholder: str,
        five_prime_flank_length: int = 8,
        three_prime_flank_length: int = 8
) -> tuple[str, str]:
    """
    Scans the FASTA file for the wildcard placeholder and dynamically extracts
    the upstream (5') and downstream (3') flanking sequences.
    """
    lhs_flank_set = set()
    rhs_flank_set = set()

    for name, seq, _ in mappy.fastx_read(fasta_path):
        idx = seq.find(bc_placeholder)
        if idx == -1:
            raise ValueError(f"Placeholder '{bc_placeholder}' not found in sequence '{name}'.")

        # Prevent negative indexing if the barcode is very close to the 5' end
        start_5p = max(0, idx - five_prime_flank_length)

        # Extract the flanks
        lhs_flank = seq[start_5p: idx]
        rhs_flank = seq[idx + len(bc_placeholder): idx + len(bc_placeholder) + three_prime_flank_length]

        lhs_flank_set.add(lhs_flank)
        rhs_flank_set.add(rhs_flank)

    if len(lhs_flank_set) != 1 or len(rhs_flank_set) != 1:
        raise ValueError(
            "Inconsistent flanks detected! Ensure all FASTA references share "
            "the exact same sequence context immediately surrounding the wildcard."
        )

    return lhs_flank_set.pop(), rhs_flank_set.pop()


# =============================================================================
# Main Pipeline
# =============================================================================

def run_pipeline(
        fastq_path: str,
        plasmid_library: dict[str, dict[str, str]],
        output_bam_path: str,
        f5_fwd: Optional[str] = None,
        f3_fwd: Optional[str] = None,
        barcode_max_edits: int = 0,
        barcode_start_pos: Optional[int] = None,
        barcode_end_pos: Optional[int] = None,
        minimap2_preset: str = "splice",
        minimap2_kmer: Optional[int] = None,
        minimap2_w: Optional[int] = None,
) -> PipelineStats:
    logger.info("Initialising aligners and barcode matcher...")

    # Automatically infer expected barcode length from the library keys
    inferred_bc_len = len(next(iter(plasmid_library)))

    # Pass the new arguments down to the builder
    aligners = build_aligner_dict(
        plasmid_library,
        preset=minimap2_preset,
        kmer=minimap2_kmer,
        w_score=minimap2_w
    )

    header, ref_map = create_bam_header(plasmid_library)

    matcher = BarcodeMatcher(sorted(plasmid_library.keys()))
    stats = PipelineStats()

    temp_bam_fd, temp_bam_path = tempfile.mkstemp(suffix=".bam")
    os.close(temp_bam_fd)

    f5_rev = mappy.revcomp(f3_fwd) if f3_fwd else None
    f3_rev = mappy.revcomp(f5_fwd) if f5_fwd else None

    try:
        logger.info("Processing reads...")
        with pysam.AlignmentFile(temp_bam_path, "wb", header=header) as unsorted_bam:
            for name, seq, qual in mappy.fastx_read(fastq_path):
                stats.processed += 1
                if stats.processed % 10_000 == 0:
                    logger.info(f"Processed {stats.processed:,} reads...")

                raw_extracted = extract_barcode_sequence(
                    seq, f5_fwd, f3_fwd, f5_rev, f3_rev,
                    expected_bc_len=inferred_bc_len,  # <-- Dynamically injected
                    start_pos=barcode_start_pos,
                    end_pos=barcode_end_pos
                )
                if not raw_extracted:
                    stats.no_flanks_found += 1
                    continue

                matched_barcode, match_result = matcher.match(
                    raw_extracted, max_edits=barcode_max_edits
                )

                if matched_barcode is None:
                    if match_result == MappingResult.BARCODE_UNRECOGNISED:
                        stats.barcode_unrecognised += 1
                    elif match_result == MappingResult.BARCODE_AMBIGUOUS:
                        stats.barcode_ambiguous += 1
                    continue

                correct_rname = plasmid_library[matched_barcode]["rname"]
                ref_id = ref_map[correct_rname]
                result = process_and_write_read(
                    unsorted_bam, ref_id, name, seq, qual, matched_barcode, aligners
                )

                if result == MappingResult.SUCCESS:
                    stats.mapped += 1
                elif result == MappingResult.MAPPING_FAILED:
                    stats.mapping_failed += 1

        logger.info("Sorting and indexing final BAM...")
        pysam.sort("-@", "4", "-o", output_bam_path, temp_bam_path)
        pysam.index("-@", "4", output_bam_path)

    finally:
        if os.path.exists(temp_bam_path):
            os.remove(temp_bam_path)

    stats.log_summary()
    return stats


# =============================================================================
# Execution
# =============================================================================

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="library-aligner: Map RNA sequencing data derived from barcoded DNA libraries.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # --- Core Inputs & Outputs ---
    core_group = parser.add_argument_group("Core Inputs & Outputs")
    core_group.add_argument("-fq", "--fastq", required=True,
                            help="Path to input FASTQ file (derived from RNA, can be gzipped).")
    core_group.add_argument("-fa", "--fasta", required=True,
                            help="Path to reference FASTA file containing template sequences.")
    core_group.add_argument("-o", "--output", required=True,
                            help="Path for the output sorted BAM file.")
    core_group.add_argument("-fq2", "--fastq2",
                            help="Path to paired read FASTQ file (for short-read sequencing).")

    # --- Barcode-to-Reference Linkage ---
    link_group = parser.add_argument_group("Barcode Definition & Linkage")
    link_group.add_argument("-c", "--csv",
                            help="Path to CSV linking barcodes to templates. Columns must be 'barcode' and 'rname'.")
    link_group.add_argument("-wc", "--wildcard",
                            help="Exact wildcard string present in the FASTA file (e.g., 'NNNNNNNNNNNNNNN').")
    link_group.add_argument("--rname_equals_barcode", action="store_true",
                            help="Flag indicating FASTA template names are exactly the barcode sequences.")

    # --- Barcode Location Details ---
    loc_group = parser.add_argument_group("Barcode Location Settings")
    loc_group.add_argument("--five_prime_flank",
                           help="Specific 5' flanking sequence (must be identical for all references).")
    loc_group.add_argument("--three_prime_flank",
                           help="Specific 3' flanking sequence (must be identical for all references).")
    loc_group.add_argument("--five_prime_flank_length", type=int, default=8,
                           help="Length of 5' flank to infer from wildcards (default: 8).")
    loc_group.add_argument("--three_prime_flank_length", type=int, default=8,
                           help="Length of 3' flank to infer from wildcards (default: 8).")
    loc_group.add_argument("--barcode_start_pos", type=int,
                           help="0-based start index of the barcode (for deterministic reads like Illumina).")
    loc_group.add_argument("--barcode_end_pos", type=int,
                           help="0-based end index of the barcode.")

    # --- Alignment & Matching Settings ---
    align_group = parser.add_argument_group("Alignment & Matching Settings")
    align_group.add_argument("--minimum_distance", type=int, default=0,
                             help="Max edit distance for fuzzy barcode matching (default: 0).")
    align_group.add_argument("--minimap2_preset", default="-ax splice",
                             help="Minimap2 preset (default: '-ax splice').")
    align_group.add_argument("--minimap2_kmer", type=int,
                             help="Minimap2 k-mer length (-k).")
    align_group.add_argument("--minimap2_w", type=int,
                             help="Minimap2 minimizer window size (-w).")

    args = parser.parse_args()

    # --- Basic Validation Logic ---
    # Ensure Option 1 requirements are met if using wildcards
    if args.wildcard and not args.csv:
        parser.error("Option 1 requires a CSV file (-c/--csv) when providing a wildcard (-wc/--wildcard).")

    # Ensure Option 2 requirements are met if using rname_equals_barcode
    if args.rname_equals_barcode:
        if args.wildcard:
            parser.error("Cannot use --wildcard and --rname_equals_barcode together.")
        if not (args.five_prime_flank and args.three_prime_flank) and not (args.barcode_start_pos is not None and args.barcode_end_pos is not None):
            parser.error("Option 2 requires explicit flanks (--five_prime_flank/--three_prime_flank) OR absolute positions (--barcode_start_pos/--barcode_end_pos).")

    try:
        # 1. Determine Flanks
        f5, f3 = None, None
        if args.five_prime_flank and args.three_prime_flank:
            print("Using explicitly provided flanking sequences...")
            f5 = args.five_prime_flank
            f3 = args.three_prime_flank
        elif args.wildcard:
            print(f"Detecting flanks using wildcard '{args.wildcard}'...")
            f5, f3 = get_flanks(
                fasta_path=args.fasta,
                bc_placeholder=args.wildcard,
                five_prime_flank_length=args.five_prime_flank_length,
                three_prime_flank_length=args.three_prime_flank_length
            )
            print(f"Detected 5' Flank: {f5}")
            print(f"Detected 3' Flank: {f3}")

        # 2. Build Library Dictionary
        print("Building library dictionary...")
        lib_dict = build_library_dictionary(
            csv_path=args.csv,
            fasta_path=args.fasta,
            to_replace=args.wildcard,
            rname_equals_barcode=args.rname_equals_barcode  # <--- Added here
        )

        # 3. Run Pipeline (You will need to pass the new args down)
        print(f"Starting pipeline. Outputting to {args.output}...")
        stats = run_pipeline(
            fastq_path=args.fastq,
            plasmid_library=lib_dict,
            output_bam_path=args.output,
            f5_fwd=f5,
            f3_fwd=f3,
            barcode_max_edits=args.minimum_distance,
            barcode_start_pos=args.barcode_start_pos,
            barcode_end_pos=args.barcode_end_pos,
            minimap2_preset=args.minimap2_preset,  # <--- Fixed
            minimap2_kmer=args.minimap2_kmer,  # <--- Fixed
            minimap2_w=args.minimap2_w  # <--- Fixed
        )
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()