"""Benchmark mappy.Aligner initialisation time.

Tests how long it takes to build an index, to see whether
on-the-fly construction (build once per read, discard after) is viable.
"""

import random
import time
import mappy

random.seed(42)

def random_seq(length: int) -> str:
    return "".join(random.choices("ACGT", k=length))


def bench_init(n_iters: int, ref_len: int, preset: str = "sr:splice", barcode_len: int = 15):
    backbone = random_seq(ref_len - barcode_len)
    insert_pos = ref_len // 2

    # Pre-generate all sequences so string ops don't pollute the timing
    seqs = []
    for _ in range(n_iters):
        bc = random_seq(barcode_len)
        seqs.append(backbone[:insert_pos] + bc + backbone[insert_pos:])

    times = []
    for seq in seqs:
        t0 = time.perf_counter()
        a = mappy.Aligner(seq=seq, preset=preset, min_chain_score=25)
        t1 = time.perf_counter()
        if not a:
            raise RuntimeError("Aligner failed")
        times.append(t1 - t0)
        del a

    avg_ms = (sum(times) / len(times)) * 1000
    min_ms = min(times) * 1000
    max_ms = max(times) * 1000
    total_s = sum(times)

    print(f"  ref_len={ref_len:>7,}  n={n_iters:>4}  |  "
          f"avg {avg_ms:>7.2f} ms  min {min_ms:>7.2f} ms  max {max_ms:>7.2f} ms  |  "
          f"total {total_s:>6.2f} s")
    return avg_ms


if __name__ == "__main__":
    print("=" * 80)
    print("mappy.Aligner initialisation benchmark")
    print("=" * 80)

    N = 200  # iterations per config

    print(f"\n--- Vary reference length (n={N} inits each, preset=splice) ---")
    for ref_len in [3_000, 5_000, 10_000, 20_000, 50_000, 100_000]:
        bench_init(N, ref_len=ref_len)

    print(f"\n--- Vary preset (ref_len=10,000, n={N}) ---")
    for preset in ["splice", "map-ont", "sr", "sr:splice"]:
        print(f"  preset={preset}")
        bench_init(N, ref_len=10_000, preset=preset)

    # Extrapolate: what would on-the-fly cost for a full run?
    print("\n--- Extrapolation: on-the-fly cost for N reads ---")
    avg = bench_init(100, ref_len=10_000)
    for n_reads in [10_000, 100_000, 1_000_000]:
        est_s = (avg / 1000) * n_reads
        est_min = est_s / 60
        print(f"    {n_reads:>10,} reads × {avg:.2f} ms  =  {est_s:>8.1f} s  ({est_min:>5.1f} min)")

    print("\n" + "=" * 80)
    print("Done.")