"""Check all twists by finite BWB computations plus a dominant-weight tail."""
import argparse
from common import Grassmannian, RHom, nonzero, require


def window(G):
    seq = [(1, -1)]
    for k in range(6):
        seq.append((0, k))
        if k <= 4:
            seq.append((1, k))
        if k <= 2:
            seq.append((2, k))
    return [(r, k, G.bundle([r, k, 0, 0, 0])) for r, k in seq]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exceptional", action="store_true",
                        help="also run the redundant exceptionality/ordering check; does not prove fullness")
    args = parser.parse_args()
    G = Grassmannian(2, 6)
    seq = window(G)
    require(len(seq) == 15, "Expected fifteen window bundles")
    pairs, cases, max_bound = 0, 0, 0
    for r, k, A in seq:
        for s, ell, B in seq:
            V = A.dual * B  # The library performs the Levi decomposition.
            weights = []
            for term in V.monomials:
                require(term.shift == 0 and term.multiplicity > 0 and len(term.factors) == 1,
                        "Expected an unshifted direct sum of irreducibles")
                d, c, *rest = term.factors[0].weight
                require(d >= 0 and not any(rest), "Expected Sym^d Udual(c)")
                weights.append((d, c))
            # Sym^d U* (c) * Sym^j U* (i) has summands
            # Sym^(d+j-2q) U* (c+i+q), 0 <= q <= min(d,j).
            # For i >= max(-c), all ambient weight coefficients are >= 0.
            # BWB then has no positive-degree cohomology, for EVERY j <= i.
            bound = max(0, *(-c for _, c in weights))
            max_bound = max(max_bound, bound)
            for i in range(bound):
                for j in range(i + 1):
                    higher = {q: n for q, n in nonzero(
                        (V * G.bundle([j, i, 0, 0, 0])).cohomology).items() if q > 0}
                    require(not higher, f"{(r,k)} -> {(s,ell)}, i={i}, j={j}: {higher}")
                    cases += 1
            pairs += 1
    print(f"PASS: {pairs} ordered pairs; {cases} finite cases; all remaining i >= j >= 0 by dominance.")
    print(f"Largest pair-specific tail threshold: i={max_bound}.")
    if args.exceptional:
        for i, (r, k, A) in enumerate(seq):
            require(nonzero(RHom(A, A)) == {0: 1}, f"Not exceptional: {(r,k)}")
            for s, ell, B in seq[:i]:
                require(not nonzero(RHom(A, B)), f"Wrong order: {(r,k)} -> {(s,ell)}")
        print("PASS: exceptionality and semiorthogonality; fullness requires the paper's separate argument.")


if __name__ == "__main__":
    main()
