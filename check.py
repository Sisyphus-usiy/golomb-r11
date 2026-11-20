"""Check exact lower bounds, complete binary covers and matching integer rulers."""
import argparse
import gzip
from pathlib import Path
import sys
from time import monotonic_ns


# These are claims that each input must prove. Only entries in 'proved' can
# supply a sub-ruler inequality. There is no import from the generator.
CLAIMS = {5: 11, 6: 17, 7: 25, 8: 34, 9: 44, 10: 55, 11: 72}


class Invalid(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise Invalid(message)


def check(stream, n, target, proved):
    line_number = 0

    def record():
        nonlocal line_number
        fields = stream.readline().split()
        line_number += 1
        require(fields, "missing record")
        return fields[0], list(map(int, fields[1:]))

    def verify_tree():
        require(record() == (b"RGC1", [n, target]), "wrong header or claim")
        tag, marks = record()
        require(tag == b"W" and len(marks) == n, "missing upper-bound witness")
        require(marks[0] == 0 and marks[-1] == target, "wrong witness endpoints")
        diffs = sorted(marks[j] - marks[i] for i in range(n) for j in range(i + 1, n))
        require(diffs[0] >= 1 and all(b - a >= 1 for a, b in zip(diffs, diffs[1:])),
                "witness fails unit separation")

        d = n - 1
        m = n * d // 2
        intervals = [(i, j) for i in range(n) for j in range(i + 1, n)]
        vectors = [tuple(int(i <= k < j) for k in range(d)) for i, j in intervals]
        # A row is its d coefficients followed by its lower bound.
        rows = [v + (m if j - i == d else proved[j - i + 1],)
                for v, (i, j) in zip(vectors, intervals)]
        rows.append((-1,) + (0,) * (d - 2) + (1, 1))  # g_last - g_first >= 1
        full_id = intervals.index((0, d))
        path = []
        nodes = leaves = empty_leaves = max_depth = 0

        def packing(tag, args):
            require(len(args) == 1 and 0 < args[0] < (1 << m), "bad interval mask")
            mask = args[0]
            ids = [i for i in range(m) if (mask >> i) & 1]
            r = len(ids)
            coeff = [sum(vectors[i][k] for i in ids) for k in range(d)]
            if tag == b"C":
                require(not (mask >> full_id) & 1, "zero included in complementary packing")
                coeff = [r - a for a in coeff]
            rows.append(tuple(coeff) + (r * (r + 1) // 2,))

        def leaf(args):
            require(len(args) >= 3 and len(args) % 2 == 1, "malformed multipliers")
            q = args[0]
            require(q >= 0, "negative length coefficient")
            total = [0] * n
            for ref, weight in zip(args[1::2], args[2::2]):
                require(weight > 0, "nonpositive multiplier")
                if ref < 0:
                    require(-ref <= len(path), "branch row outside the current path")
                    row = path[-1 - ref]
                else:
                    require(ref < len(rows), "undeclared row")
                    row = rows[ref]
                for k, a in enumerate(row):
                    total[k] += weight * a
            require(total[:d] == [q] * d, "coefficient identity failed")
            require(total[d] >= q * target if q else total[d] > 0, "insufficient bound")
            return q == 0

        def node():
            nonlocal nodes, leaves, empty_leaves, max_depth
            tag, args = record()
            while tag in (b"P", b"C"):
                packing(tag, args)
                tag, args = record()
            nodes += 1
            max_depth = max(max_depth, len(path))
            if tag == b"L":
                empty_leaves += leaf(args)
                leaves += 1
                return
            require(tag == b"B" and len(args) == 2, "expected a split or a leaf")
            u, v = args
            require(0 <= u < m and 0 <= v < m and u != v, "invalid difference comparison")
            delta = tuple(a - b for a, b in zip(vectors[u], vectors[v]))
            path.append(delta + (1,))
            node()
            path[-1] = tuple(-a for a in delta) + (1,)
            node()
            path.pop()

        node()
        require(stream.read(1) == b"", "trailing data after the root tree")
        require(nodes == 2 * leaves - 1, "incomplete binary tree")
        return nodes, leaves, empty_leaves, max_depth

    try:
        return verify_tree()
    except (ValueError, RecursionError) as exc:
        raise Invalid(f"n={n}, record {line_number}: {exc}") from exc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", type=Path, default=Path(__file__).resolve().parent / "certificates")
    parser.add_argument("--through", type=int, choices=range(5, 12), default=11)
    args = parser.parse_args()
    start = monotonic_ns()
    # Sorted unit-separated differences have order statistics at least 1,2,...
    proved = {n: n * (n - 1) // 2 for n in (2, 3, 4)}
    for n in range(5, args.through + 1):
        step = monotonic_ns()
        with gzip.open(args.dir / f"n{n}.rgc.gz", "rb") as stream:
            nodes, leaves, empty, depth = check(stream, n, CLAIMS[n], proved)
        proved[n] = CLAIMS[n]
        ms = (monotonic_ns() - step) // 1_000_000
        print(f"PASS n={n} G_R={proved[n]} nodes={nodes} leaves={leaves} "
              f"empty={empty} max_depth={depth} check_ms={ms}", flush=True)
    print(f"PASS through n={args.through} total_ms={(monotonic_ns() - start) // 1_000_000}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, EOFError) as exc:
        sys.exit(f"REJECT: {exc}")
