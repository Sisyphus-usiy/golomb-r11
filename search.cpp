// Generate a binary chamber cover and exact LP witnesses. Build as C++17.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using Int = std::int64_t;
using Mask = std::uint64_t;
constexpr Int state_limit = 10'000'000;
constexpr std::array<int, 12> targets{0, 0, 1, 3, 6, 11, 17, 25, 34, 44, 55, 72};
const std::array<std::vector<int>, 12> witnesses{{
    {}, {}, {0,1}, {0,1,3}, {0,1,4,6}, {0,1,4,9,11}, {0,1,4,10,12,17},
    {0,1,4,10,18,23,25}, {0,1,4,9,15,22,32,34}, {0,1,5,12,25,27,35,41,44},
    {0,1,6,10,23,26,34,41,53,55}, {0,1,4,13,28,33,47,54,64,70,72}
}};

int triangular(int r) { return r * (r + 1) / 2; }

// Bounds on every stored basis entry make all intermediates fit in int64.
// For d <= 10 we have |a_j| <= 55 and b <= 1540. In a pivot the largest
// possible multiply-add has magnitude <= 2640 * state_limit^2 < 2^63.
Int checked(Int value) {
    if (value < -state_limit || value > state_limit)
        throw std::overflow_error("basis exceeds the checked arithmetic range");
    return value;
}

Int quotient(Int numerator, Int denominator) {
    if (denominator == 1) return checked(numerator);
    if (numerator % denominator != 0) throw std::logic_error("inexact pivot");
    return checked(numerator / denominator);
}

template<int D>
class Search {
    static constexpr int N = D + 1;
    static constexpr int M = N * D / 2;
    using Vec = std::array<Int, D>;
    using Coeff = std::array<int, D>;
    struct Row { Coeff a{}; int rhs = 0; int ref = 0; };
    struct Interval { int lo, hi; };
    struct Pair { int u, v; };
    struct Basis {
        // B is the matrix of tight rows. inv = den * B^-1, den = |det B|.
        Int den = 1;
        std::array<Vec, D> inv{};
        Vec x{}, dual{};  // Scaled primal point and dual multipliers.
        std::array<int, D> row{};
    };
    struct Cut { Mask mask = 0; bool complement = false; Int violation = 0; };

    std::ostream& out;
    const int target;
    std::array<Interval, M> intervals{};
    std::vector<Row> active;
    std::vector<Pair> comparisons;
    std::unordered_map<Mask, int> packing_ids;
    Int nodes = 0, leaves = 0, empty = 0, pivots = 0;
    int max_depth = 0;

    static Int dot(const Coeff& a, const Vec& x) {
        return std::inner_product(a.begin(), a.end(), x.begin(), Int{0});
    }

    Row order_row(int u, int v, int ref) const {
        Row row;
        row.rhs = 1;
        row.ref = ref;
        for (int j = 0; j < D; ++j) row.a[j] = active[u].a[j] - active[v].a[j];
        return row;
    }

    Basis root() {
        Basis b;
        for (int lo = 0; lo < D; ++lo) {
            for (int hi = lo + 1; hi <= D; ++hi) {
                int id = static_cast<int>(active.size());
                intervals[id] = {lo, hi};
                Row row;
                row.ref = id;
                std::fill(row.a.begin() + lo, row.a.begin() + hi, 1);
                row.rhs = hi - lo == D ? M : targets[hi - lo + 1];
                active.push_back(row);
                if (hi == lo + 1) b.row[lo] = id;
            }
            b.inv[lo][lo] = b.x[lo] = b.dual[lo] = 1;
        }
        active.push_back(order_row(M - 1, 0, M));  // last gap >= first gap + 1

        // Overlapping intervals cancel to two disjoint intervals; containment
        // already separates their lengths. Short supports win exact ties.
        for (int u = 0; u < M; ++u)
            for (int v = u + 1; v < M; ++v)
                if (intervals[u].hi <= intervals[v].lo) comparisons.push_back({u, v});
        auto key = [&](Pair p) {
            auto a = intervals[p.u], c = intervals[p.v];
            return std::make_tuple(a.hi - a.lo + c.hi - c.lo, c.hi - a.lo, p.u, p.v);
        };
        std::sort(comparisons.begin(), comparisons.end(),
                  [&](Pair a, Pair c) { return key(a) < key(c); });
        return b;
    }

    void length_leaf(const Basis& b) {
        Int gcd = b.den;
        for (Int y : b.dual) gcd = std::gcd(gcd, y);
        out << "L " << b.den / gcd;
        for (int i = 0; i < D; ++i) {
            if (b.dual[i] < 0) throw std::logic_error("dual feasibility lost");
            if (b.dual[i]) out << ' ' << active[b.row[i]].ref << ' ' << b.dual[i] / gcd;
        }
        out << '\n';
        ++leaves;
    }

    void empty_leaf(const Basis& b, int entering, const Vec& transform) {
        out << "L 0 " << active[entering].ref << ' ' << b.den;
        for (int i = 0; i < D; ++i)
            if (transform[i]) out << ' ' << active[b.row[i]].ref << ' ' << -transform[i];
        out << '\n';
        ++empty;
        ++leaves;
    }

    // Enter a violated row in the ordinary simplex for the dual
    // max b^t y, A^t y = 1, y >= 0. The minimum ratio preserves y >= 0.
    bool pivot(Basis& b, int entering, Int deficiency) {
        const auto& a = active[entering].a;
        Vec t{};
        for (int j = 0; j < D; ++j)
            for (int i = 0; i < D; ++i) t[i] += a[j] * b.inv[j][i];
        int k = -1;
        for (int i = 0; i < D; ++i) {
            if (t[i] <= 0) continue;
            if (k < 0 || b.dual[i] * t[k] < b.dual[k] * t[i] ||
                (b.dual[i] * t[k] == b.dual[k] * t[i] && b.row[i] < b.row[k])) k = i;
        }
        if (k < 0) {
            empty_leaf(b, entering, t);
            return false;
        }
        const Int h = checked(t[k]);
        for (int j = 0; j < D; ++j) {
            const Int column = b.inv[j][k];
            b.x[j] = quotient(h * b.x[j] + deficiency * column, b.den);
            for (int i = 0; i < D; ++i)
                if (i != k) b.inv[j][i] = quotient(h * b.inv[j][i] - column * t[i], b.den);
        }
        for (int i = 0; i < D; ++i)
            if (i != k) b.dual[i] = quotient(h * b.dual[i] - b.dual[k] * t[i], b.den);
        b.den = h;
        b.row[k] = entering;
        ++pivots;
        return true;
    }

    // Return true if the node closed. A dual bound can close it before a
    // primal feasible optimum is needed. Bland's rule handles degeneracy.
    bool optimize(Basis& b) {
        for (int iteration = 0; ; ++iteration) {
            if (std::accumulate(b.x.begin(), b.x.end(), Int{0}) >= target * b.den) {
                length_leaf(b);
                return true;
            }
            Int best = 0;
            int entering = -1;
            for (int i = 0; i < static_cast<int>(active.size()); ++i) {
                const Row& r = active[i];
                Int deficiency = r.rhs * b.den - dot(r.a, b.x);
                if (deficiency > best) {
                    best = deficiency;
                    entering = i;
                    if (iteration >= 24) break;
                }
            }
            if (entering < 0) return false;
            if (!pivot(b, entering, best)) return true;
        }
    }

    std::array<Int, M> differences(const Basis& b) const {
        std::array<Int, N> marks{};
        std::partial_sum(b.x.begin(), b.x.end(), marks.begin() + 1);
        std::array<Int, M> values{};
        for (int i = 0; i < M; ++i)
            values[i] = marks[intervals[i].hi] - marks[intervals[i].lo];
        return values;
    }

    Cut separate(const std::array<Int, M>& values, Int den) const {
        std::array<int, M> order{};
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int i, int j) {
            return std::tie(values[i], i) < std::tie(values[j], j);
        });
        Cut best;
        // The full ruler is the unique longest interval, so omit it from
        // the list of complementary distances measured down from the length.
        for (bool complement : {false, true}) {
            Int sum = 0;
            Mask mask = 0;
            for (int r = 1; r <= M - int(complement); ++r) {
                int id = complement ? order[M - 1 - r] : order[r - 1];
                sum += complement ? values[order.back()] - values[id] : values[id];
                mask |= Mask{1} << id;
                Int violation = triangular(r) * den - sum;
                if (violation > best.violation) best = {mask, complement, violation};
            }
        }
        return best;
    }

    void add_cut(const Cut& cut) {
        Mask key = cut.mask | (Mask(cut.complement) << 63);
        auto result = packing_ids.emplace(key, M + 1 + static_cast<int>(packing_ids.size()));
        if (result.second) out << (cut.complement ? "C " : "P ") << cut.mask << '\n';
        Row row;
        row.ref = result.first->second;
        int count = 0;
        for (int i = 0; i < M; ++i) {
            if (!(cut.mask & (Mask{1} << i))) continue;
            ++count;
            for (int j = 0; j < D; ++j)
                row.a[j] += cut.complement ? 1 - active[i].a[j] : active[i].a[j];
        }
        row.rhs = triangular(count);
        active.push_back(row);
    }

    void visit(Basis b, int depth) {
        ++nodes;
        max_depth = std::max(max_depth, depth);
        if (nodes % 1'000'000 == 0) std::cerr << "visited " << nodes << " nodes\n";
        const auto saved = active.size();
        while (!optimize(b)) {
            auto values = differences(b);
            Cut cut = separate(values, b.den);
            if (cut.violation > 0) {
                add_cut(cut);
                continue;
            }
            Pair split{-1, -1};
            Int best = b.den;
            for (Pair p : comparisons) {
                Int distance = std::abs(values[p.u] - values[p.v]);
                if (distance < best) {
                    best = distance;
                    split = p;
                    if (distance == 0) break;
                }
            }
            if (split.u < 0) {
                Int mark = 0;
                std::cerr << "shorter ruler: 0";
                for (Int gap : b.x) {
                    mark += gap;
                    std::cerr << ' ' << mark << '/' << b.den;
                }
                std::cerr << '\n';
                throw std::runtime_error("open chamber contains a shorter ruler");
            }
            out << "B " << split.u << ' ' << split.v << '\n';
            active.push_back(order_row(split.u, split.v, -1 - depth));
            visit(b, depth + 1);
            active.back() = order_row(split.v, split.u, -1 - depth);
            visit(b, depth + 1);
            break;
        }
        active.resize(saved);
    }

public:
    Search(std::ostream& stream, int requested) : out(stream), target(requested) {}

    void run() {
        out << "RGC1 " << N << ' ' << target << "\nW";
        for (int mark : witnesses[N]) out << ' ' << mark;
        out << '\n';
        auto b = root();
        visit(b, 0);
        if (nodes != 2 * leaves - 1) throw std::logic_error("tree count mismatch");
        std::cerr << "n=" << N << " nodes=" << nodes << " leaves=" << leaves
                  << " empty=" << empty << " cuts=" << packing_ids.size()
                  << " max_depth=" << max_depth << " pivots=" << pivots << '\n';
    }
};

int main(int argc, char** argv) {
    try {
        if (argc != 3 && argc != 4) throw std::runtime_error("usage: search N OUTPUT [TARGET]");
        int n = std::stoi(argv[1]);
        if (n < 5 || n > 11) throw std::runtime_error("N must be 5 to 11");
        int target = argc == 4 ? std::stoi(argv[3]) : targets[n];
        if (target < 1 || target > 100) throw std::runtime_error("TARGET must be 1 to 100");
        std::ofstream out(argv[2], std::ios::binary);
        out.exceptions(std::ios::badbit | std::ios::failbit);
        auto start = std::chrono::steady_clock::now();
        switch (n) {
        case 5: Search<4>(out, target).run(); break;
        case 6: Search<5>(out, target).run(); break;
        case 7: Search<6>(out, target).run(); break;
        case 8: Search<7>(out, target).run(); break;
        case 9: Search<8>(out, target).run(); break;
        case 10: Search<9>(out, target).run(); break;
        case 11: Search<10>(out, target).run(); break;
        }
        out.close();
        std::cerr << "search_ms=" << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start).count() << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
