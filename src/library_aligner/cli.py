import argparse
import logging
import sys
from .core import get_flanks, build_library_dictionary, run_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Fast, chimera-safe Nanopore barcode extraction and plasmid alignment."
    )
    parser.add_argument("-q", "--fastq", required=True, help="Path to input FASTQ file (can be gzipped)")
    parser.add_argument("-r", "--refs", required=True, help="Path to reference multi-FASTA")
    parser.add_argument("-c", "--csv", required=True, help="Path to barcode-reference linkage CSV")
    parser.add_argument("-o", "--out", required=True, help="Path for output sorted BAM file")
    parser.add_argument("-p", "--placeholder", default="NNNNNNNNNNNNNNN",
                        help="Barcode placeholder in FASTA (default: 15 Ns)")
    parser.add_argument("-e", "--max-edits", type=int, default=0,
                        help="Max edit distance for barcode matching (default: 0)")

    args = parser.parse_args()

    try:
        print("Detecting flanking sequences...")
        f5, f3 = get_flanks(args.refs, args.placeholder, flank_len=10)

        print("Building library dictionary...")
        lib_dict = build_library_dictionary(args.csv, args.refs, to_replace=args.placeholder)

        print(f"Starting pipeline. Outputting to {args.out}...")
        stats = run_pipeline(
            fastq_path=args.fastq,
            plasmid_library=lib_dict,
            output_bam_path=args.out,
            f5_fwd=f5,
            f3_fwd=f3,
            barcode_max_edits=args.max_edits
        )
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()