"""Profile memory usage of mappy.Aligner instances.

Generates synthetic plasmid references of varying sizes and counts,
then measures RSS growth per aligner to estimate real-world costs.
"""

import os
import random
import psutil
import mappy

random.seed(42)

def random_seq(length: int) -> str:
    return "".join(random.choices("ACGT", k=length))

def rss_mb() -> float:
    return psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)

def profile(n_aligners: int, ref_len: int, preset: str = "splice", barcode_len: int = 15):
    """Build n_aligners, each with a unique reference of ref_len bp (barcode region randomised)."""

    # Generate a shared backbone with a barcode slot
    backbone = random_seq(ref_len - barcode_len)
    insert_pos = ref_len // 2

    before = rss_mb()
    aligners = []
    for i in range(n_aligners):
        barcode = random_seq(barcode_len)
        seq = backbone[:insert_pos] + barcode + backbone[insert_pos:]
        a = mappy.Aligner(seq=seq, preset=preset, min_chain_score=25)
        if not a:
            raise RuntimeError(f"Aligner {i} failed")
        aligners.append(a)

    after = rss_mb()
    total = after - before
    per_aligner = total / n_aligners

    print(f"  {n_aligners:>6,} aligners × {ref_len:>7,} bp  |  "
          f"total {total:>8.1f} MB  |  per-aligner {per_aligner:>6.2f} MB")

    # Keep reference alive until we're done measuring
    del aligners
    return total, per_aligner


def profile_shared_backbone(n_barcodes: int, ref_len: int, preset: str = "splice"):
    """Measure cost of the single-aligner alternative: one index for the shared backbone."""
    backbone = random_seq(ref_len)

    before = rss_mb()
    a = mappy.Aligner(seq=backbone, preset=preset, min_chain_score=25)
    if not a:
        raise RuntimeError("Shared aligner failed")
    after = rss_mb()

    total = after - before
    print(f"  {'1 (shared)':>14} × {ref_len:>7,} bp  |  "
          f"total {total:>8.1f} MB  |  (serves {n_barcodes:,} barcodes)")
    del a
    return total


if __name__ == "__main__":
    print("=" * 72)
    print("mappy.Aligner memory profile")
    print("=" * 72)

    # --- Vary library size with fixed reference length ---
    print("\n--- Scaling by library size (ref_len = 10,000 bp) ---")
    for n in [10, 50, 100, 250, 500]:
        profile(n, ref_len=10_000)

    # --- Vary reference length with fixed library size ---
    print("\n--- Scaling by reference length (n = 100 aligners) ---")
    for ref_len in [3_000, 10_000, 30_000]:
        profile(100, ref_len=ref_len)

    # --- Compare: shared backbone alternative ---
    print("\n--- Shared backbone alternative ---")
    for n in [100, 500]:
        profile(n, ref_len=10_000)
        profile_shared_backbone(n, ref_len=10_000)
        print()

    print("=" * 72)
    print("Done.")