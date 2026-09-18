"""Exact degree support on U/Gamma: finite pushforwards and a proved tail."""
from functools import lru_cache
from common import Grassmannian, nonzero, require

P = Grassmannian(1, 6)


@lru_cache(maxsize=None)
def base_degrees(i, m):
    # Partition (m; i,0,0,0,0) corresponds to [m-i,i,0,0,0].
    return frozenset(nonzero(P.bundle([m-i, i, 0, 0, 0]).cohomology))


@lru_cache(maxsize=None)
def degrees(a, b):
    """Return exact nonzero degrees, not dimensions of the infinite direct sum.

    Relative O(t) on P(O(-h)+3O) pushes to Sym^t(O(h)+3O)
    for t>=0; for t<=-4 it pushes in degree 3 to
    Sym^(-t-4)(O(-h)+3O) * O(-h). Intermediate twists vanish.
    Positive multiplicities do not affect degree support.
    """
    cutoff = max(0, -b, -a-1, -a-b)
    result = set()
    for i in range(cutoff):
        t = i + b
        if t >= 0:
            for j in range(t + 1):
                result.update(base_degrees(i, i+a+j))
        elif t <= -4:
            for j in range(-t-3):
                result.update(q+3 for q in base_degrees(i, i+a-j-1))
    # For i>=cutoff, m=i+a+j >= -1. The shifted GL6 weight is
    # (m+5,i+4,3,2,1,0): degree 0 when m>=i; degree 1 when
    # -1<=m<=i-2; singular when m=i-1. No other degrees occur.
    # j=max(0,-a) witnesses degree 0, and j=0 witnesses degree 1
    # precisely when a<=-2. The cutoff makes both witnesses admissible.
    require(base_degrees(cutoff, cutoff+a+max(0, -a)) == {0}, "H0 tail witness failed")
    result.add(0)
    if a <= -2:
        require(base_degrees(cutoff, cutoff+a) == {1}, "H1 tail witness failed")
        result.add(1)
    return frozenset(result)


CLAIMS = [
    (0, -5, 5, {0}), (1, -5, 6, {0, 3}), (1, -3, 6, {0}),
    (-1, -6, 5, {0, 4}), (-1, -4, 5, {0}),
    (2, -4, 6, {0, 1, 3}), (2, -3, 6, {0, 1}),
    (-2, -6, 4, {0, 1, 4}), (-2, -5, 4, {0, 1}),
    (3, -2, 5, {0}), (-3, -5, 2, {0, 1}),
    (4, 0, 4, {0}), (-4, -4, 0, {0, 1}),
]


def main():
    print(" a / b " + "".join(f"{b:>8}" for b in range(-6, 7)))
    for a in range(-4, 5):
        print(f"{a:>6} " + "".join(f"{','.join(map(str, sorted(degrees(a,b)))):>8}"
                                     for b in range(-6, 7)))
    for a, lo, hi, allowed in CLAIMS:
        for b in range(lo, hi+1):
            require(degrees(a, b) <= allowed, f"Corollary fails at {(a,b)}: {degrees(a,b)}")
    print("PASS: all 13 ranges in the nine corollary items; no arbitrary truncation.")
    print("H0 is infinite dimensional for all (a,b); H1 is infinite dimensional exactly for a <= -2.")


if __name__ == "__main__":
    main()
