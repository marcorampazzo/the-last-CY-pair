"""Reproduce both appendix Hodge diamonds using the current library API."""
from common import Grassmannian, require
from hodge import hodge_numbers, _zero_locus_omega_bounds, _tighten_hodge_bounds
from varieties import ProjectiveBundle, ZeroLocus


def blowup_hodge_numbers(Z):
    bounds, chi = {}, {}
    for p in range(Z.dimension + 1):
        row, chi[p] = _zero_locus_omega_bounds(Z, p)
        bounds.update({(p, q): interval for q, interval in row.items()})
    lo, hi = bounds[0, 0]
    bounds[0, 0] = (max(1, lo), hi)
    _tighten_hodge_bounds(bounds, chi, Z.dimension)
    print("Bounds before the blow-up input:",
          {pq: interval for pq, interval in bounds.items() if interval[0] != interval[1]})
    # Z is the blow-up of a smooth projective
    # threefold along a nonempty smooth curve. The exceptional divisor adds
    # an independent (1,1) class to the pullback of an ample class: h11>=2.
    # This uses neither Y's Calabi-Yau property nor the target Hodge numbers.
    print("Geometric input from prop:gorCY3: a blow-up along a nonempty smooth curve gives h11 >= 2.")
    lo, hi = bounds[1, 1]
    require(hi >= 2, "The computed upper bound contradicts the blow-up input")
    bounds[1, 1] = (max(2, lo), hi)
    _tighten_hodge_bounds(bounds, chi, Z.dimension)
    require(all(lo == hi for lo, hi in bounds.values()), f"Unresolved Hodge bounds: {bounds}")
    h = {pq: lo for pq, (lo, _) in bounds.items()}
    for p in range(4):
        require(sum((-1 if q % 2 else 1) * h[p, q] for q in range(4)) == chi[p],
                f"Euler row {p} failed")
        for q in range(4):
            require(h[p, q] >= 0 and h[p, q] == h[q, p] == h[3-p, 3-q],
                    "Hodge symmetry or Serre duality failed")
    print("Exact chi(Omega^p):", chi)
    return h


def main():
    G = Grassmannian(2, 6)
    X = ZeroLocus(G, G.Udual * G.O(1) + 3 * G.O(1))
    P = Grassmannian(1, 6)
    W = ProjectiveBundle(P, P.O(-1) + 3 * P.trivial_bundle)
    Yt = ZeroLocus(W, P.Qdual * P.O(1) * W.O(1))
    for name, Z, expected in [("X", X, (1, 52)), ("Ytilde", Yt, (2, 53))]:
        print(f"Computing {name} (smooth regular section assumed)...", flush=True)
        require(Z.dimension == 3, f"{name}: unexpected dimension")
        h = hodge_numbers(Z) if name == "X" else blowup_hodge_numbers(Z)
        require((h[1, 1], h[1, 2]) == expected, f"{name}: unexpected Hodge numbers: {h}")
        require(h[0, 0] == h[3, 0] == 1 and h.get((1, 0), 0) == h.get((2, 0), 0) == 0,
                f"{name}: unexpected structure-sheaf cohomology")
        euler = sum((-1 if (p + q) % 2 else 1) * v for (p, q), v in h.items())
        print(f"dim={Z.dimension}; canonical bundle={Z.canonical_bundle}")
        print(f"h11={h[1, 1]}; h12={h[1, 2]}; topological Euler characteristic={euler}")
        for degree in range(7):
            print(" ".join(str(h.get((p, degree-p), 0)) for p in range(4)
                           if 0 <= degree-p <= 3).center(24))


if __name__ == "__main__":
    main()
