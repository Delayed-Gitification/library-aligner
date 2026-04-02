"""library-aligner: Map RNA sequencing data derived from barcoded DNA libraries.

Extracts barcodes from reads (via flanking sequences or fixed positions),
matches them to a known library, then aligns each read against the
barcode-specific reference and writes a sorted, indexed BAM file.
"""

import os
import sys
import tempfile
import logging
from dataclasses import dataclass
from enum import Enum
from typing import Optional
import mappy
import edlib
import pysam
import pandas as pd
import argparse
import pybktree

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
    def __init__(self, known_barcodes: list[str], check_reverse_complement: bool = True, fuzzy_mode: str = "NW"):
        self.known_barcodes = known_barcodes
        self.exact_set = set(known_barcodes)
        self.check_reverse_complement = check_reverse_complement

        def _dist(a, b):
            return edlib.align(a, b, mode=fuzzy_mode, task="distance")["editDistance"]

        # Build once at init — O(N log N), paid once
        self.tree = pybktree.BKTree(_dist, known_barcodes)

        # Always initialise RC structures to safe defaults
        self.rc_to_original = {}
        self.rc_tree = None

        if check_reverse_complement:
            self.rc_to_original = {mappy.revcomp(bc): bc for bc in known_barcodes}
            if len(self.rc_to_original) != len(known_barcodes):
                raise ValueError("Two barcodes in the library are reverse complements of each other.")
            self.rc_tree = pybktree.BKTree(_dist, list(self.rc_to_original.keys()))

    def match(self, extracted_seq: str, max_edits: int = 2) -> tuple[Optional[str], MappingResult]:
        rc_extracted = mappy.revcomp(extracted_seq) if self.check_reverse_complement else None

        # 1. Fast path: exact match — O(1)
        if extracted_seq in self.exact_set:
            return extracted_seq, MappingResult.SUCCESS
        if self.check_reverse_complement and rc_extracted in self.exact_set:
            return rc_extracted, MappingResult.SUCCESS

        if max_edits < 1:
            return None, MappingResult.BARCODE_UNRECOGNISED

        # 2. Fuzzy path: BK-tree query — O(log N) average
        fwd_hits = self.tree.find(extracted_seq, max_edits)  # [(dist, bc), ...]
        rev_hits = self.rc_tree.find(rc_extracted, max_edits) if self.check_reverse_complement else []

        all_hits = list(fwd_hits)  # [(dist, bc) for dist, bc in fwd_hits]
        all_hits += [(dist, self.rc_to_original[bc]) for dist, bc in rev_hits]

        if not all_hits:
            return None, MappingResult.BARCODE_UNRECOGNISED

        min_dist = min(dist for dist, _ in all_hits)
        best = {bc for dist, bc in all_hits if dist == min_dist}

        if len(best) > 1:
            return None, MappingResult.BARCODE_AMBIGUOUS

        return best.pop(), MappingResult.SUCCESS


def extract_barcode_sequence(
        read_seq: str, f5_fwd=None, f3_fwd=None,
        expected_bc_len=15, max_error_rate=0.2,
        start_pos=None, end_pos=None,
        search_first_n_bases=None, search_last_n_bases=None,
        length_tolerance=2,
        check_reverse_complement: bool = True,
) -> list[str]:
    candidates = []

    if not f5_fwd and not f3_fwd:
        if start_pos is not None and end_pos is not None and end_pos <= len(read_seq):
            candidates.append(read_seq[start_pos:end_pos])
        return candidates

    mode = "dual"
    if f5_fwd and not f3_fwd:
        mode = "f5_only"
    elif f3_fwd and not f5_fwd:
        mode = "f3_only"

    if mode == "dual":
        fwd_template = f5_fwd + ("N" * expected_bc_len) + f3_fwd
        max_errs = int((len(f5_fwd) + len(f3_fwd)) * max_error_rate) + length_tolerance
    elif mode == "f5_only":
        fwd_template = f5_fwd + ("N" * expected_bc_len)
        max_errs = int(len(f5_fwd) * max_error_rate)
    elif mode == "f3_only":
        fwd_template = ("N" * expected_bc_len) + f3_fwd
        max_errs = int(len(f3_fwd) * max_error_rate)

    dna_eq = [("N", "A"), ("N", "C"), ("N", "G"), ("N", "T")]

    read_states = [("f", read_seq)]
    if check_reverse_complement:
        read_states.append(("r", mappy.revcomp(read_seq)))

    for orient, current_seq in read_states:
        windows_to_check = []

        if search_first_n_bases is not None:
            windows_to_check.append(current_seq[:search_first_n_bases])
        if search_last_n_bases is not None:
            windows_to_check.append(current_seq[-search_last_n_bases:])

        if not windows_to_check:
            windows_to_check.append(current_seq)

        for w_idx, window_seq in enumerate(windows_to_check):

            # THE SCAFFOLD SEARCH
            res = edlib.align(fwd_template, window_seq, "HW", "locations", max_errs, additionalEqualities=dna_eq)

            if res["editDistance"] == -1:
                continue

            s, e = res["locations"][0]
            footprint = window_seq[max(0, s - 3): min(len(window_seq), e + 4)]

            if mode == "dual":
                f5_errs = int(len(f5_fwd) * max_error_rate)
                f3_errs = int(len(f3_fwd) * max_error_rate)
                res_5 = edlib.align(f5_fwd, footprint, "HW", "locations", k=f5_errs)
                res_3 = edlib.align(f3_fwd, footprint, "HW", "locations", k=f3_errs)

                if res_5["editDistance"] != -1 and res_3["editDistance"] != -1:
                    e5, s3 = res_5["locations"][0][1], res_3["locations"][-1][0]
                    if s3 > e5:
                        extracted = footprint[e5 + 1: s3]
                        if abs(len(extracted) - expected_bc_len) <= length_tolerance:
                            candidates.append(extracted)

            elif mode == "f5_only":
                f5_errs = int(len(f5_fwd) * max_error_rate)
                res_5 = edlib.align(f5_fwd, footprint, "HW", "locations", k=f5_errs)

                if res_5["editDistance"] != -1:
                    e5 = res_5["locations"][0][1]
                    extracted = footprint[e5 + 1: e5 + 1 + expected_bc_len]

                    if len(extracted) == expected_bc_len:
                        candidates.append(extracted)

            elif mode == "f3_only":
                f3_errs = int(len(f3_fwd) * max_error_rate)
                res_3 = edlib.align(f3_fwd, footprint, "HW", "locations", k=f3_errs)

                if res_3["editDistance"] != -1:
                    s3 = res_3["locations"][-1][0]
                    extracted = footprint[max(0, s3 - expected_bc_len): s3]

                    if len(extracted) == expected_bc_len:
                        candidates.append(extracted)

    return candidates


# =============================================================================
# BAM Handling
# =============================================================================


def process_and_write_read(bam_writer, ref_id, name, seq1, qual1, barcode, aligners_dict,
                           seq2=None, qual2=None, min_mapq=0):
    """Map reads and write BAM records. min_mapq filters alignments below this
    mapping quality threshold (0 = accept all)."""

    aligner = aligners_dict.get(barcode)
    if not aligner:
        raise RuntimeError(
            f"No aligner found for barcode '{barcode}' — this is a bug, barcode should have been validated before reaching this point.")

    alns1 = [a for a in aligner.map(seq1) if a.mapq >= min_mapq]
    if not alns1:
        return MappingResult.MAPPING_FAILED

    q_arr1 = pysam.qualitystring_to_array(qual1) if qual1 else None

    # Pysam CIGAR ops: M=0, I=1, D=2, N=3, S=4, H=5, P=6, =7, X=8
    _CIGAR_OPS = "MIDNSHP=X"

    def build_cigar(aln, seq_len):
        # mappy yields (length, op); pysam expects (op, length)
        cigar = [(op, length) for length, op in aln.cigar]
        clip_l, clip_r = aln.q_st, seq_len - aln.q_en
        if aln.strand == 1:
            if clip_l > 0: cigar.insert(0, (4, clip_l))
            if clip_r > 0: cigar.append((4, clip_r))
        else:
            if clip_r > 0: cigar.insert(0, (4, clip_r))
            if clip_l > 0: cigar.append((4, clip_l))
        return cigar

    def cigar_to_string(cigar_tuples):
        """Convert pysam-style (op, length) tuples to a SAM CIGAR string."""
        return "".join(f"{length}{_CIGAR_OPS[op]}" for op, length in cigar_tuples)

    def build_sa_tag(alns, seq_len, ref_name):
        """Build SA:Z tag value for a list of alignments.
        Format per entry: rname,pos,strand,CIGAR,mapQ,NM;
        """
        entries = []
        for a in alns:
            strand = "-" if a.strand == -1 else "+"
            cigar_str = cigar_to_string(build_cigar(a, seq_len))
            # NM tag: use mappy's NM if available, else 0
            nm = getattr(a, "NM", 0)
            entries.append(f"{ref_name},{a.r_st + 1},{strand},{cigar_str},{a.mapq},{nm}")
        return ";".join(entries) + ";"

    def set_sequence(seg, seq, q_arr):
        """Set query sequence and qualities on a pysam AlignedSegment.
        Must be called BEFORE setting cigar."""
        if seg.is_reverse:
            seg.query_sequence = mappy.revcomp(seq)
            if q_arr is not None:
                seg.query_qualities = q_arr[::-1]
        else:
            seg.query_sequence = seq
            if q_arr is not None:
                seg.query_qualities = q_arr

    # Resolve reference name for SA tags
    ref_name = bam_writer.header.references[ref_id]

    # --- Paired-end path ---
    if seq2 is not None:
        alns2 = [a for a in aligner.map(seq2) if a.mapq >= min_mapq]
        q_arr2 = pysam.qualitystring_to_array(qual2) if qual2 else None

        # Best primary alignment from each read — used for mate coordinate fields
        aln1 = alns1[0]
        aln2 = alns2[0] if alns2 else None

        # SAM spec TLEN: rightmost mapped base of downstream mate minus
        # leftmost mapped base of upstream mate, +1, with sign indicating orientation.
        if aln2 is not None:
            leftmost = min(aln1.r_st, aln2.r_st)
            rightmost = max(aln1.r_en, aln2.r_en)
            tlen_abs = rightmost - leftmost
            # Positive for the leftmost mate, negative for the rightmost
            tlen = tlen_abs if aln1.r_st <= aln2.r_st else -tlen_abs
        else:
            tlen = 0

        for aln in alns1:
            a = pysam.AlignedSegment(bam_writer.header)
            a.query_name = name
            a.reference_id = ref_id
            a.reference_start = aln.r_st
            a.mapping_quality = aln.mapq
            a.is_paired = True
            a.is_read1 = True
            a.is_reverse = (aln.strand == -1)
            a.is_supplementary = not aln.is_primary

            if aln2 is not None:
                a.next_reference_id = ref_id
                a.next_reference_start = aln2.r_st
                a.mate_is_reverse = (aln2.strand == -1)
                a.is_proper_pair = aln.is_primary  # only flag proper pair on primary
                a.template_length = tlen
            else:
                a.next_reference_id = ref_id
                a.next_reference_start = aln.r_st
                a.mate_is_unmapped = True

            # Sequence → qualities → cigar
            set_sequence(a, seq1, q_arr1)
            a.cigar = build_cigar(aln, len(seq1))
            # SA tag: list other alignments for this read (supplementary support)
            if len(alns1) > 1:
                other_alns = [x for x in alns1 if x is not aln]
                a.set_tag("SA", build_sa_tag(other_alns, len(seq1), ref_name))
            bam_writer.write(a)

        for aln in (alns2 if alns2 else []):
            a = pysam.AlignedSegment(bam_writer.header)
            a.query_name = name
            a.reference_id = ref_id
            a.reference_start = aln.r_st
            a.mapping_quality = aln.mapq
            a.is_paired = True
            a.is_read2 = True
            a.is_reverse = (aln.strand == -1)
            a.is_supplementary = not aln.is_primary
            a.next_reference_id = ref_id
            a.next_reference_start = aln1.r_st
            a.mate_is_reverse = (aln1.strand == -1)
            a.is_proper_pair = aln.is_primary
            a.template_length = -tlen  # R2 gets the negated R1 TLEN

            set_sequence(a, seq2, q_arr2)
            a.cigar = build_cigar(aln, len(seq2))
            if len(alns2) > 1:
                other_alns = [x for x in alns2 if x is not aln]
                a.set_tag("SA", build_sa_tag(other_alns, len(seq2), ref_name))
            bam_writer.write(a)

        # Unmapped R2 placeholder
        if not alns2:
            a = pysam.AlignedSegment(bam_writer.header)
            a.query_name = name
            a.is_paired = True
            a.is_read2 = True
            a.is_unmapped = True
            a.reference_id = ref_id
            a.reference_start = aln1.r_st  # convention: use mate's position
            a.next_reference_id = ref_id
            a.next_reference_start = aln1.r_st
            a.mate_is_reverse = (aln1.strand == -1)
            a.query_sequence = seq2
            if q_arr2 is not None:
                a.query_qualities = q_arr2
            bam_writer.write(a)

        return MappingResult.SUCCESS

    # --- Single-end path ---
    for aln in alns1:
        a = pysam.AlignedSegment(bam_writer.header)
        a.query_name = name
        a.reference_id = ref_id
        a.reference_start = aln.r_st
        a.mapping_quality = aln.mapq
        a.is_reverse = (aln.strand == -1)
        a.is_supplementary = not aln.is_primary

        set_sequence(a, seq1, q_arr1)
        a.cigar = build_cigar(aln, len(seq1))
        if len(alns1) > 1:
            other_alns = [x for x in alns1 if x is not aln]
            a.set_tag("SA", build_sa_tag(other_alns, len(seq1), ref_name))
        bam_writer.write(a)

    return MappingResult.SUCCESS

# =============================================================================
# Setup Helpers
# =============================================================================

class LazyAlignerDict:
    """Drop-in replacement for the aligner dict that builds each aligner on
    demand and discards it after use.  Trades CPU for memory: ~0.8 ms per
    build for a 10 kb reference, but holds only one aligner at a time."""

    def __init__(self, plasmid_references, preset, kmer=None, w_score=None):
        self._refs = plasmid_references
        self._preset = preset
        self._kmer = kmer
        self._w = w_score

    def get(self, barcode):
        ref_data = self._refs.get(barcode)
        if ref_data is None:
            return None
        kwargs = {"seq": ref_data["seq"], "preset": self._preset, "min_chain_score": 25}
        if self._kmer is not None: kwargs["k"] = self._kmer
        if self._w is not None: kwargs["w"] = self._w
        aligner = mappy.Aligner(**kwargs)
        if not aligner:
            raise RuntimeError(f"Failed to build on-the-fly aligner for barcode '{barcode}'")
        return aligner


def build_aligner_dict(plasmid_references, preset="splice", kmer=None, w_score=None, cache=True):
    clean_preset = preset.replace("-ax ", "").replace("-x ", "").strip()

    if not cache:
        logger.info("Aligner caching disabled (--no_mappy_cache). Aligners will be built on the fly.")
        return LazyAlignerDict(plasmid_references, clean_preset, kmer, w_score)

    aligners = {}
    seq_to_aligner = {}  # keyed on actual sequence — safe regardless of wildcard usage

    for barcode, ref_data in plasmid_references.items():
        seq = ref_data["seq"]
        if seq not in seq_to_aligner:
            kwargs = {"seq": seq, "preset": clean_preset, "min_chain_score": 25}
            if kmer is not None: kwargs["k"] = kmer
            if w_score is not None: kwargs["w"] = w_score
            aligner = mappy.Aligner(**kwargs)
            if not aligner:
                raise RuntimeError(f"Failed to load reference for barcode '{barcode}' using preset '{clean_preset}'")
            seq_to_aligner[seq] = aligner
        aligners[barcode] = seq_to_aligner[seq]

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
) -> tuple[Optional[str], Optional[str]]:
    lhs_flank_set = set()
    rhs_flank_set = set()

    for name, seq, _ in mappy.fastx_read(fasta_path):
        idx = seq.find(bc_placeholder)
        if idx == -1:
            raise ValueError(f"Placeholder '{bc_placeholder}' not found in sequence '{name}'.")

        start_5p = max(0, idx - five_prime_flank_length)

        lhs_flank = seq[start_5p: idx]
        rhs_flank = seq[idx + len(bc_placeholder): idx + len(bc_placeholder) + three_prime_flank_length]

        lhs_flank_set.add(lhs_flank)
        rhs_flank_set.add(rhs_flank)

    if len(lhs_flank_set) != 1 or len(rhs_flank_set) != 1:
        raise ValueError(
            "Inconsistent flanks detected! Ensure all FASTA references share "
            "the exact same sequence context immediately surrounding the wildcard."
        )

    lhs = lhs_flank_set.pop()
    rhs = rhs_flank_set.pop()

    # Return None instead of empty strings if the barcode is at the absolute edge
    return (lhs if lhs else None), (rhs if rhs else None)

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
        check_reverse_complement: bool = True,
        fastq2_path: Optional[str] = None,       # <-- new
        barcode_read: int = 1,
        barcode_length_tolerance: int = 2,
        fuzzy_mode: str = 'HW',
        cache_aligners: bool = True,
        search_first_n_bases: Optional[int] = None,  # <-- Added
        search_last_n_bases: Optional[int] = None
) -> PipelineStats:
    logger.info("Initialising aligners and barcode matcher...")

    # Automatically infer expected barcode length from the library keys
    inferred_bc_len = len(next(iter(plasmid_library)))

    aligners = build_aligner_dict(
        plasmid_library,
        preset=minimap2_preset,
        kmer=minimap2_kmer,
        w_score=minimap2_w,
        cache=cache_aligners,
    )

    header, ref_map = create_bam_header(plasmid_library)

    # Pass the flag into the matcher
    matcher = BarcodeMatcher(
        sorted(plasmid_library.keys()),
        check_reverse_complement=check_reverse_complement,
        fuzzy_mode=fuzzy_mode
    )
    stats = PipelineStats()

    temp_bam_fd, temp_bam_path = tempfile.mkstemp(suffix=".bam")
    os.close(temp_bam_fd)

    f5_rev = mappy.revcomp(f3_fwd) if f3_fwd else None
    f3_rev = mappy.revcomp(f5_fwd) if f5_fwd else None

    try:
        logger.info("Processing reads...")
        with pysam.AlignmentFile(temp_bam_path, "wb", header=header) as unsorted_bam:

            read_iter = mappy.fastx_read(fastq_path)
            read2_iter = mappy.fastx_read(fastq2_path) if fastq2_path else None

            for r1 in read_iter:
                name1, seq1, qual1 = r1

                if read2_iter is not None:
                    r2 = next(read2_iter, None)
                    if r2 is None:
                        logger.warning("R2 exhausted before R1 — FASTQs may be mismatched.")
                        break
                    name2, seq2, qual2 = r2
                    # Strip /1 /2 or .1 .2 suffixes for comparison
                    base1 = name1.rsplit("/", 1)[0].rsplit(" ", 1)[0]
                    base2 = name2.rsplit("/", 1)[0].rsplit(" ", 1)[0]
                    if base1 != base2:
                        logger.error(
                            f"Read name mismatch at read {stats.processed + 1}: "
                            f"R1='{name1}' vs R2='{name2}'. FASTQs are out of sync."
                        )
                        raise RuntimeError("Paired FASTQ files are out of sync — read names do not match.")
                else:
                    seq2, qual2 = None, None

                stats.processed += 1
                if stats.processed % 10_000 == 0:
                    logger.info(f"Processed {stats.processed:,} reads...")

                # Route barcode extraction to the correct read
                bc_seq = seq1 if (barcode_read == 1 or seq2 is None) else seq2

                # Notice this is now explicitly expecting a LIST
                raw_extracted_list = extract_barcode_sequence(
                    bc_seq, f5_fwd, f3_fwd,
                    expected_bc_len=inferred_bc_len,
                    start_pos=barcode_start_pos,
                    end_pos=barcode_end_pos,
                    search_first_n_bases=search_first_n_bases,
                    search_last_n_bases=search_last_n_bases,
                    length_tolerance=barcode_length_tolerance,
                    check_reverse_complement=check_reverse_complement,
                )

                if not raw_extracted_list:
                    stats.no_flanks_found += 1
                    continue

                matched_barcode = None
                match_result = MappingResult.NO_FLANKS_FOUND

                # Test every candidate extracted from the read
                for raw_extracted in raw_extracted_list:
                    bc, res = matcher.match(raw_extracted, max_edits=barcode_max_edits)
                    if bc is not None:
                        matched_barcode = bc
                        match_result = res
                        break  # Stop looking once we find a valid barcode!
                    else:
                        match_result = res  # Store the last failure reason for stats

                if matched_barcode is None:
                    if match_result == MappingResult.BARCODE_UNRECOGNISED:
                        stats.barcode_unrecognised += 1
                    elif match_result == MappingResult.BARCODE_AMBIGUOUS:
                        stats.barcode_ambiguous += 1
                    continue

                # --- Proceed with alignment as normal ---
                correct_rname = plasmid_library[matched_barcode]["rname"]
                ref_id = ref_map[correct_rname]
                result = process_and_write_read(
                    unsorted_bam, ref_id, name1, seq1, qual1, matched_barcode, aligners,
                    seq2=seq2, qual2=qual2
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
    loc_group.add_argument("--barcode_read", type=int, choices=[1, 2], default=1,
                           help="Which read the barcode is on when using positional extraction (default: 1).")
    loc_group.add_argument("--barcode_length_tolerance", type=int, default=2,
                           help="Max allowed deviation in bp from expected barcode length during flank-based extraction (default: 2).")
    loc_group.add_argument("--search_first_n_bases", type=int,
                           help="Number of bases to scan at the 5' end of the read (e.g., 350).")
    loc_group.add_argument("--search_last_n_bases", type=int,
                           help="Number of bases to scan at the 3' end of the read (e.g., 350).")

    # --- Alignment & Matching Settings ---
    align_group = parser.add_argument_group("Alignment & Matching Settings")
    align_group.add_argument("--check_reverse_complement", choices=["True", "False"], default="True",
                             help="Check the reverse complement of the barcode (default: True). Must be False for positional extraction.")
    align_group.add_argument("--minimum_distance", type=int, default=0,
                             help="Max edit distance for fuzzy barcode matching (default: 0).")
    align_group.add_argument("--minimap2_preset", default="-ax splice",
                             help="Minimap2 preset (default: '-ax splice').")
    align_group.add_argument("--minimap2_kmer", type=int,
                             help="Minimap2 k-mer length (-k).", default=10)
    align_group.add_argument("--minimap2_w", type=int,
                             help="Minimap2 minimizer window size (-w).", default=4)
    align_group.add_argument("--fuzzy_mode", choices=["NW", "HW"], default=None,
                             help="Edit distance mode for fuzzy barcode matching. NW=global (Illumina), HW=infix (ONT). Auto-inferred if not set.")
    align_group.add_argument("--no_mappy_cache", action="store_true",
                             help="Build minimap2 aligners on the fly instead of caching all in memory. Slower but uses minimal RAM.")

    args = parser.parse_args()

    check_rc = args.check_reverse_complement == "True"

    if args.fuzzy_mode is not None:
        fuzzy_mode = args.fuzzy_mode
    else:
        # Positional extraction = fixed-length barcode = global distance is correct
        # Flank-based extraction = variable-length extraction = infix distance is safer
        fuzzy_mode = "NW" if args.barcode_start_pos is not None else "HW"

    # --- Basic Validation Logic ---
        # Check if the user has provided AT LEAST one flank
        has_any_flank = bool(args.five_prime_flank or args.three_prime_flank)

        is_exact_positional = not args.wildcard and not has_any_flank and (args.barcode_start_pos is not None)

        if args.wildcard and not args.csv:
            parser.error("Option 1 requires a CSV file (-c/--csv) when providing a wildcard (-wc/--wildcard).")

        if is_exact_positional and check_rc:
            parser.error(
                "When using exact positional extraction (no flanks/wildcards), you must set --check_reverse_complement False.")

        if args.barcode_read == 2 and not args.fastq2:
            parser.error("--barcode_read 2 requires a paired FASTQ (--fastq2).")

        if args.rname_equals_barcode:
            if args.wildcard:
                parser.error("Cannot use --wildcard and --rname_equals_barcode together.")
            if not has_any_flank and not is_exact_positional:
                parser.error(
                    "Option 2 requires at least one flank (--five_prime_flank and/or --three_prime_flank) OR exact absolute positions (--barcode_start_pos/--barcode_end_pos).")

    try:
        # 1. Determine Flanks
        # 1. Determine Flanks
        f5 = args.five_prime_flank
        f3 = args.three_prime_flank

        # If the user didn't provide ANY explicit flanks, but did provide a wildcard, auto-detect them
        if not (f5 or f3) and args.wildcard:
            print(f"Detecting flanks using wildcard '{args.wildcard}'...")
            f5, f3 = get_flanks(
                fasta_path=args.fasta,
                bc_placeholder=args.wildcard,
                five_prime_flank_length=args.five_prime_flank_length,
                three_prime_flank_length=args.three_prime_flank_length
            )
            print(f"Detected 5' Flank: {f5}")
            print(f"Detected 3' Flank: {f3}")
        else:
            print("Using explicitly provided flanking sequences (Single or Dual)...")

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
            minimap2_preset=args.minimap2_preset,
            minimap2_kmer=args.minimap2_kmer,
            minimap2_w=args.minimap2_w,
            check_reverse_complement=check_rc,
            fastq2_path=args.fastq2,  # <-- new
            barcode_read=args.barcode_read,  # <-- new
            barcode_length_tolerance=args.barcode_length_tolerance,
            fuzzy_mode=fuzzy_mode,
            cache_aligners=not args.no_mappy_cache,
            search_first_n_bases=args.search_first_n_bases,  # <-- Added
            search_last_n_bases=args.search_last_n_bases,
        )
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
