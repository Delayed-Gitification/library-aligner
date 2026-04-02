import time
import random
import edlib
import mappy
import pybktree
from typing import Optional


def generate_distinct_barcodes(n: int, length: int = 25, min_edit_dist: int = 4) -> list[str]:
    bases = ['A', 'C', 'G', 'T']
    accepted = []
    attempts = 0
    while len(accepted) < n:
        attempts += 1
        if attempts > 100_000:
            raise RuntimeError(f"Could not generate {n} barcodes after 100k attempts.")
        candidate = "".join(random.choices(bases, k=length))
        if all(
            True
            # edlib.align(candidate, bc, mode="NW", task="distance")["editDistance"] >= min_edit_dist
            for bc in accepted
        ):
            accepted.append(candidate)
    return accepted


def add_noise(seq: str, error_rate: float) -> str:
    bases = ['A', 'C', 'G', 'T']
    return "".join(
        random.choice([b for b in bases if b != base]) if random.random() < error_rate else base
        for base in seq
    )


# --- Linear matcher (old approach) ---
class LinearBarcodeMatcher:
    def __init__(self, known_barcodes):
        self.known_barcodes = known_barcodes
        self.exact_set = set(known_barcodes)

    def match(self, extracted_seq: str, max_edits: int = 1) -> Optional[str]:
        if extracted_seq in self.exact_set:
            return extracted_seq
        best_matches = set()
        min_dist = max_edits + 1
        for bc in self.known_barcodes:
            res = edlib.align(bc, extracted_seq, mode="NW", task="distance", k=max_edits)
            if res["editDistance"] != -1:
                dist = res["editDistance"]
                if dist < min_dist:
                    min_dist, best_matches = dist, {bc}
                elif dist == min_dist:
                    best_matches.add(bc)
        return best_matches.pop() if len(best_matches) == 1 else None


# --- BK-tree matcher (new approach) ---
class BKTreeBarcodeMatcher:
    def __init__(self, known_barcodes):
        self.exact_set = set(known_barcodes)

        def _dist(a, b):
            return edlib.align(a, b, mode="NW", task="distance")["editDistance"]

        self.tree = pybktree.BKTree(_dist, known_barcodes)

    def match(self, extracted_seq: str, max_edits: int = 1) -> Optional[str]:
        if extracted_seq in self.exact_set:
            return extracted_seq
        hits = self.tree.find(extracted_seq, max_edits)
        if not hits:
            return None
        min_dist = min(d for d, _ in hits)
        best = {bc for d, bc in hits if d == min_dist}
        return best.pop() if len(best) == 1 else None


def test_matcher_benchmark():
    random.seed(42)

    library_sizes = [50, 200, 500, 1000, 10_000]
    n_queries = 2000       # queries per benchmark run
    noise_rate = 0.05      # high enough that most queries miss exact match and hit fuzzy path

    print("\n\n" + "=" * 80)
    print("BARCODE MATCHER BENCHMARK: Linear O(N) vs BK-Tree O(log N)")
    print("=" * 80)
    print(f"Queries per run: {n_queries:,} | Noise rate: {noise_rate:.0%} | Max edits: 1")
    print(f"(Exact matches excluded from timing — both are O(1) hash lookup)\n")

    print(f"  {'Library size':<16} {'Linear (ms)':<18} {'BK-Tree (ms)':<18} {'Speedup':<12} {'Results match'}")
    print(f"  {'-'*76}")

    for n in library_sizes:
        barcodes = generate_distinct_barcodes(n, length=25, min_edit_dist=1)

        # Generate queries that will miss exact match and exercise the fuzzy path
        queries = [add_noise(random.choice(barcodes), noise_rate) for _ in range(n_queries)]
        # Remove any accidental exact matches so we're only timing the fuzzy path
        queries = [q for q in queries if q not in set(barcodes)]

        linear  = LinearBarcodeMatcher(barcodes)
        bktree  = BKTreeBarcodeMatcher(barcodes)

        # Warm up
        for q in queries[:10]:
            linear.match(q)
            bktree.match(q)

        # Time linear
        t0 = time.perf_counter()
        linear_results = [linear.match(q) for q in queries]
        t_linear = (time.perf_counter() - t0) * 1000

        # Time BK-tree
        t0 = time.perf_counter()
        bktree_results = [bktree.match(q) for q in queries]
        t_bktree = (time.perf_counter() - t0) * 1000

        speedup = t_linear / t_bktree
        results_match = linear_results == bktree_results

        print(f"  {n:<16} {t_linear:<18.1f} {t_bktree:<18.1f} {speedup:<12.1f} {'✅' if results_match else '❌ MISMATCH'}")

    print("=" * 80 + "\n")