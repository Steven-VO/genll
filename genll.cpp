// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Steven Van Overberghe
//
// genll: generates the locally linear graphs on N vertices (every edge in exactly one
// triangle), one per isomorphism class.  See README.md.
//
// Compiled with -DMAXN=WORDSIZE -DWORDSIZE=32 (genll, N <= 32) or 64 (genllL, N <= 64)
// against nautyW1.a or nautyL1.a.
//
// Canonical construction path:
//   augmentation  add a vertex u (the hub) and k >= 0 disjoint pairs {v_i,w_i} of old
//                 vertices, with edges u v_i, u w_i, v_i w_i.  Feasible when
//                   (P1) d(v_i,w_i) >= 3,
//                   (P2) no edge between two chosen vertices.
//                 C4-free: d(v_i,w_i) >= 4, and d(x,y) >= 3 for x,y in different pairs.
//   reduction     delete a vertex and the matching on its neighbourhood.
//   canonical     the vertex of minimum degree; ties broken by the invariant kappa and
//                 then by a canonical labelling from nauty.
//   orbits        only the lexicographically least pair set of each Aut(G)-orbit is
//                 expanded.
//
// -5 rejects children with a 5-cycle through a new edge.  C5-free locally linear
// graphs are also C4-free.
//

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "nauty.h"
#include "nautyBitset.h"
#include "graph.h"

static_assert(MAXN == WORDSIZE, "link against the MAXN=WORDSIZE nauty archives");

constexpr int MaxOrder = WORDSIZE;

typedef InverseBitset<WORDSIZE> Set;
typedef BitGraph<Set, MaxOrder> Graph;

// The adjacency rows are passed to nauty as they are.
static_assert(sizeof(Set) == sizeof(setword));

typedef std::array<int, MaxOrder> Perm;

// A pair set a1 b1 a2 b2 ... (a_i < b_i, a_1 < a_2 < ...) packed into one integer, so
// that comparing integers compares pair sets lexicographically.
typedef unsigned __int128 PairKey;
constexpr int PairBits = (MaxOrder <= 32) ? 5 : 6;
constexpr int MaxPairsInKey = 128 / (2 * PairBits);

// kappa(v) = (degree, near, far, farther, colour), where B_r(v) is the ball of radius r:
//   near     |B_2(v)|, sum of deg(y)^2 over N(v)
//   far      |B_3(v)|, sum of deg(y) over B_2(v)
//   farther  |B_4(v)|, sum of deg(y) over B_3(v)
//   colour   colour refinement starting from the degrees
// In tryChild, each component is only computed for the vertices still tied with the hub.
constexpr int SquareBits = 18;         // sum of deg^2 <= 64 * 63^2
constexpr int BallBits   = 7;          // the size of a ball is at most 64
constexpr int DegreeSumBits = 13;      // sum of deg over a ball <= 64 * 63
constexpr int NearBits   = BallBits + SquareBits;

// rounds of colour refinement in kappa
#ifndef REFINEMENT_ROUNDS
#define REFINEMENT_ROUNDS 3
#endif

// statistics (-DGENLL_STATS=1); the counters are empty when disabled
#ifndef GENLL_STATS
#define GENLL_STATS 0
#endif

template<bool On>
struct CounterOf {
    struct Absent {};
    [[no_unique_address]] std::conditional_t<On, long long, Absent> value {};
    void operator++(int) { if constexpr (On) ++value; }
    void raiseTo(size_t v) {
        if constexpr (On) if (static_cast<long long>(v) > value) value = static_cast<long long>(v);
    }
    long long get() const { if constexpr (On) return value; else return -1; }
};
typedef CounterOf<GENLL_STATS != 0> Counter;

// res/mod
struct SplitState {
    long long counter = 0;
    int residue = 0, modulus = 1, level = 0;
};

// -d#: per level, the vertices that must be chosen and the least number of pairs
struct MinDegreeState {
    int minimum = 0;
    Set forcedAt[MaxOrder + 1] {};
    int leastPairsAt[MaxOrder + 1] {};
};

// -M: two vertices are far if a new triangle could contain both (distance >= 3, or >= 4
// for C4-free).  A graph is maximal iff no three vertices are pairwise far.  On the last
// level, every far triangle of the parent must be destroyed by the augmentation; a
// greedy family of disjoint far triangles is used to prune the pair search.
constexpr int MaxFarTriangles = MaxOrder / 3 + 1;

struct MaximalState {
    Set triangle[MaxFarTriangles] {};          // for bounding the search
    uint8_t corner[MaxFarTriangles][3] {};     // the same three vertices, for testing
    int triangles = 0;
};

struct Absent {};

static inline uint64_t mixColour(uint64_t x) {
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}

// nauty callback: collect the generators of the automorphism group
static std::vector<Perm>* automorphismTarget = nullptr;

static void collectAutomorphism(int, int* perm, int*, int, int, int n) {
    automorphismTarget->resize(automorphismTarget->size() + 1);
    Perm& image = automorphismTarget->back();
    for (int i = 0; i < n; ++i) image[i] = perm[i];
}

// ---------------------------------------------------------------------------

template<bool C4Free, bool Split, bool CountOnly, bool MinDegree, bool Maximal>
struct Generator {
    int targetOrder = 0;
    bool c5free = false;           // -5

    // vertices at distance > FarRadius are far
    static constexpr int FarRadius = C4Free ? 3 : 2;

    [[no_unique_address]] std::conditional_t<Split, SplitState, Absent> split {};
    [[no_unique_address]] std::conditional_t<MinDegree, MinDegreeState, Absent> minDegree {};
    [[no_unique_address]] std::conditional_t<Maximal, MaximalState, Absent> maximal {};

    // --- the graph under construction -------------------------------------
    Graph g {};                    // vertices 0..order-1, plus the hub while building
    int degree[MaxOrder] {};

    // The arrays below are indexed by level, i.e. by the order of the parent.

    // B_2(v) and the sum of squared neighbour degrees, in the parent
    std::array<Graph, MaxOrder + 1> closed2 {};
    uint32_t squaresAt[MaxOrder + 1][MaxOrder] {};

    // partner of each chosen vertex in the current pair set
    int partnerAt[MaxOrder + 1][MaxOrder] {};

    // B_3(v) in the parent (C4-free only)
    std::conditional_t<C4Free, std::array<Graph, MaxOrder + 1>, std::array<Graph, 0>> closed3 {};

    // the pair set being built
    uint8_t pairsAt[MaxOrder + 1][MaxOrder] {};

    // maximum number of pairs, the vertices a pair set of that size must contain, and
    // the vertices of each degree
    int maxChooseAt[MaxOrder + 1] {};
    Set mustCoverAt[MaxOrder + 1] {};
    Set degreeClassAt[MaxOrder + 1][MaxOrder + 3] {};

    // generators of Aut(G)
    std::vector<Perm> generators[MaxOrder + 1];

    // --- nauty ------------------------------------------------------------
    DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    int lab[MaxOrder], ptn[MaxOrder], orbits[MaxOrder];
    ::graph canonical[MaxOrder];

    // --- the orbit of one pair set ----------------------------------------
    std::vector<PairKey> orbit;

    // --- output ------------------------------------------------------------
    static constexpr int OutputBufferSize = 1 << 16;
    static constexpr int LongestGraph6 = longestGraph6(MaxOrder);
    char outputBuffer[OutputBufferSize];
    int outputUsed = 0;

    // --- statistics -------------------------------------------------------
    long long outCount = 0;
    Counter nodeCount, childCount, nautyCount, canonCount, orbitTestCount, largestOrbit;
    Counter degreeRejects, orbitRejects, kappaRejects;
    Counter shortcutAccepts, nautyAccepts, nautyRejects;
    Counter maximalRejects, triangleRejects;

    Generator() {
        options.defaultptn = FALSE;
        for (auto& generator : generators) generator.reserve(2 * MaxOrder);
        orbit.reserve(1024);
    }

    // =======================================================================
    // The pair set as an integer
    // =======================================================================

    static PairKey packPairs(const uint8_t* sequence, int k) {
        PairKey key = 0;
        for (int i = 0; i < 2 * k; ++i) key = (key << PairBits) | sequence[i];
        return key;
    }

    static void unpackPairs(PairKey key, int k, uint8_t* sequence) {
        for (int i = 2 * k - 1; i >= 0; --i) {
            sequence[i] = uint8_t(key & ((PairKey{1} << PairBits) - 1));
            key >>= PairBits;
        }
    }

    // sort each pair, then the pairs by their first vertex
    static void normalisePairs(uint8_t* sequence, int k) {
        for (int i = 0; i < k; ++i)
            if (sequence[2 * i] > sequence[2 * i + 1])
                std::swap(sequence[2 * i], sequence[2 * i + 1]);
        for (int i = 1; i < k; ++i) {
            uint8_t a = sequence[2 * i], b = sequence[2 * i + 1];
            int j = i - 1;
            while (j >= 0 && sequence[2 * j] > a) {
                sequence[2 * j + 2] = sequence[2 * j];
                sequence[2 * j + 3] = sequence[2 * j + 1];
                --j;
            }
            sequence[2 * j + 2] = a;
            sequence[2 * j + 3] = b;
        }
    }

    // Is the current pair set the least of its Aut(G)-orbit?  The orbit is enumerated
    // breadth first and the search stops at the first smaller image.
    bool isLeastInOrbit(int order, int k) {
        const std::vector<Perm>& gens = generators[order];
        if (k == 0 || gens.empty()) return true;

        orbitTestCount++;
        const PairKey start = packPairs(pairsAt[order], k);
        orbit.clear();
        orbit.push_back(start);

        uint8_t sequence[MaxOrder];
        for (size_t head = 0; head < orbit.size(); ++head) {
            for (const Perm& perm : gens) {
                unpackPairs(orbit[head], k, sequence);
                for (int i = 0; i < 2 * k; ++i) sequence[i] = uint8_t(perm[sequence[i]]);
                normalisePairs(sequence, k);
                const PairKey image = packPairs(sequence, k);
                if (image < start) return false;
                if (std::find(orbit.begin(), orbit.end(), image) == orbit.end())
                    orbit.push_back(image);
            }
        }
        largestOrbit.raiseTo(orbit.size());
        return true;
    }

    // =======================================================================
    // The canonical reduction
    // =======================================================================

    // The near component of kappa in the child, also returning B_2.  Computed from the
    // parent's closed2 and squaresAt: only the hub and the chosen vertices changed.
    uint64_t hubNearInvariant(int order, Set hubBall, uint64_t hubSquares, Set& ball) const {
        ball = hubBall | Set::bitAt(order);
        const uint64_t result = (uint64_t(ball.popcount()) << SquareBits) | hubSquares;
#ifdef GENLL_VERIFY_KAPPA
        verifyNear(order, order, result, ball);
#endif
        return result;
    }

    uint64_t nearInvariant(int order, int x, Set& ball) const {
        const Set chosen = g[order];
        ball = closed2[order][x];
        uint64_t squares = squaresAt[order][x];

        const int* partner = partnerAt[order];
        if (chosen.hasBit(x)) {
            const int p = partner[x];
            ball |= chosen | Set::bitAt(order)
                  | (g[p] & ~(Set::bitAt(order) | Set::bitAt(x)));
            squares += uint64_t(degree[order]) * degree[order]
                     + uint64_t(degree[p]) * degree[p];
        } else {
            const Set touched = g[x] & chosen;
            if (touched.isNotEmpty()) {
                ball |= Set::bitAt(order);
                for (int y : touched) {
                    ball |= Set::bitAt(partner[y]);
                    squares += 4 * uint64_t(degree[y]) - 4;
                }
            }
        }
        const uint64_t result = (uint64_t(ball.popcount()) << SquareBits) | squares;
#ifdef GENLL_VERIFY_KAPPA
        verifyNear(order, x, result, ball);
#endif
        return result;
    }

    // For a ball B_r: |B_{r+1}| and the degree sum over B_r; B_{r+1} is returned.
    uint64_t ballInvariant(Set ball, Set& wider) const {
        wider = ball;
        uint64_t degreeSum = 0;
        for (int y : ball) { wider |= g[y]; degreeSum += degree[y]; }
        return (uint64_t(wider.popcount()) << DegreeSumBits) | degreeSum;
    }

    // direct computation, for -DGENLL_VERIFY_KAPPA
    uint64_t nearInvariantDirect(int v, Set& ball) const {
        const Set neighbours = g[v];
        uint64_t squares = 0;
        ball = neighbours | Set::bitAt(v);
        for (int y : neighbours) {
            squares += uint64_t(degree[y]) * degree[y];
            ball |= g[y];
        }
        return (uint64_t(ball.popcount()) << SquareBits) | squares;
    }

    void verifyNear(int order, int v, uint64_t got, Set gotBall) const {
        Set wantBall;
        const uint64_t want = nearInvariantDirect(v, wantBall);
        if (want != got || !(wantBall == gotBall)) {
            std::cerr << "genll: incremental kappa disagrees at order " << order
                      << " vertex " << v << "\n";
            exit(2);
        }
    }

    // The colour component of kappa: colour refinement with hashed colours.  A hash
    // collision only makes the invariant weaker.
    void refineColours(int order, uint64_t* colour) const {
        for (int x = 0; x < order; ++x) colour[x] = degree[x];
        uint64_t next[MaxOrder];
        for (int round = 0; round < REFINEMENT_ROUNDS; ++round) {
            for (int x = 0; x < order; ++x) {
                uint64_t sum = 0;
                for (int y : g[x]) sum += mixColour(colour[y]);
                next[x] = mixColour(colour[x]) * 0x9E3779B97F4A7C15ull + sum;
            }
            for (int x = 0; x < order; ++x) colour[x] = next[x];
        }
    }

    // Run nauty on the child, for the canonical labelling and/or the generators.  The
    // first cell of the colouring must be exactly the kappa-minimal vertices (the hub
    // and 'tied'); the other cells only use degree and near.
    void runNauty(int order, bool wantCanonical, bool wantGenerators,
                  Set hubBall, uint64_t hubSquares, const int* tied, int ties) {
        const int hub = order - 1;
        uint64_t nearKey[MaxOrder];
        Set ball;
        for (int x = 0; x < order; ++x) {
            const uint64_t near = (x == hub)
                    ? hubNearInvariant(hub, hubBall, hubSquares, ball)
                    : nearInvariant(hub, x, ball);
            nearKey[x] = (uint64_t(degree[x]) << NearBits) | near;
            lab[x] = x;
        }
        for (int i = 1; i < order; ++i) {
            const int v = lab[i];
            const uint64_t k = nearKey[v];
            int j = i - 1;
            while (j >= 0 && nearKey[lab[j]] > k) { lab[j + 1] = lab[j]; --j; }
            lab[j + 1] = v;
        }
        for (int i = 0; i + 1 < order; ++i) ptn[i] = (nearKey[lab[i]] == nearKey[lab[i + 1]]);
        ptn[order - 1] = 0;

        // split the kappa-minimal vertices off the first cell
        int first = 0;
        while (ptn[first]) ++first;
        ++first;                                   // the first near class is lab[0..first)
        if (ties + 1 < first) {
            Set minimal = Set::bitAt(hub);
            for (int i = 0; i < ties; ++i) minimal.add(tied[i]);
            int front = 0, back = ties + 1, sorted[MaxOrder];
            for (int i = 0; i < first; ++i)
                (minimal.hasBit(lab[i]) ? sorted[front++] : sorted[back++]) = lab[i];
            for (int i = 0; i < first; ++i) lab[i] = sorted[i];
            ptn[ties] = 0;
        }

        options.getcanon = wantCanonical;
        options.userautomproc = wantGenerators ? collectAutomorphism : nullptr;
        if (wantGenerators) {
            generators[order].clear();
            automorphismTarget = &generators[order];
        }

        nautyCount++;
        if (wantCanonical) canonCount++;
        densenauty(reinterpret_cast<::graph*>(g.adj.data()), lab, ptn, orbits,
                   &options, &stats, 1, order, wantCanonical ? canonical : nullptr);
    }

    // =======================================================================
    // One candidate child
    // =======================================================================

    // -5: is there a path x-a-b-c-y on five distinct vertices?
    bool pathOfLengthFour(int x, int y) const {
        const Set notXY = ~(Set::bitAt(x) | Set::bitAt(y));
        const Set ys = g[y] & ~Set::bitAt(x);
        for (int a : g[x] & ~Set::bitAt(y)) {
            for (int c : ys & ~Set::bitAt(a))
                if ((g[a] & g[c] & notXY).isNotEmpty()) return true;
        }
        return false;
    }

    // -5: is there a 5-cycle through a new edge?
    bool hasNewC5(int hub) const {
        for (int v : g[hub]) {
            if (pathOfLengthFour(hub, v)) return true;
            const int w = partnerAt[hub][v];
            if (v < w && pathOfLengthFour(v, w)) return true;
        }
        return false;
    }

    void tryChild(int order, int k, Set hubBall, uint64_t hubSquares) {
        childCount++;
        const int hub = order;
        if (c5free && hasNewC5(hub)) return;
        const int childOrder = order + 1;
        const int hubDegree = 2 * k;

        // The hub needs minimum degree.  This can only fail for the maximum k, where all
        // vertices of minimum degree must be chosen.
        if (k == maxChooseAt[order] && (mustCoverAt[order] & ~g[hub]).isNotEmpty()) {
            degreeRejects++;
            return;
        }
        // -M, last level: every far triangle must be destroyed
        if constexpr (Maximal) {
            if (childOrder == targetOrder) {
                for (int i = 0; i < maximal.triangles; ++i) {
                    const uint8_t* corner = maximal.corner[i];
                    if (!closeInChild(corner[0], corner[1])
                        && !closeInChild(corner[0], corner[2])
                        && !closeInChild(corner[1], corner[2])) { triangleRejects++; return; }
                }
            }
        }

        if (!isLeastInOrbit(order, k)) { orbitRejects++; return; }

        // rivals: the vertices with the same degree as the hub in the child
        const Set chosen = g[hub];
        Set rivals = degreeClassAt[order][hubDegree] & ~chosen;
        if (hubDegree >= 2) rivals |= degreeClassAt[order][hubDegree - 2] & chosen;

        // The hub must have the least kappa among them.  If no rival ties with it, no
        // canonical labelling is needed.
        Set hubBallWhole, ball;
        int tied[MaxOrder], ties = 0;
        if (rivals.isNotEmpty()) {
            const uint64_t hubNear = hubNearInvariant(hub, hubBall, hubSquares, hubBallWhole);
            for (int x : rivals) {
                const uint64_t near = nearInvariant(order, x, ball);
                if (near < hubNear) { kappaRejects++; return; }
                if (near == hubNear) tied[ties++] = x;
            }
        }
        Set tiedBall[MaxOrder], hubWider, wider;
        if (ties > 0) {
            const uint64_t hubFar = ballInvariant(hubBallWhole, hubWider);
            int remaining = 0;
            for (int i = 0; i < ties; ++i) {
                nearInvariant(order, tied[i], ball);
                const uint64_t far = ballInvariant(ball, wider);
                if (far < hubFar) { kappaRejects++; return; }
                if (far == hubFar) { tiedBall[remaining] = wider; tied[remaining++] = tied[i]; }
            }
            ties = remaining;
        }
        if (ties > 0) {
            const uint64_t hubFarther = ballInvariant(hubWider, wider);
            int remaining = 0;
            for (int i = 0; i < ties; ++i) {
                const uint64_t farther = ballInvariant(tiedBall[i], wider);
                if (farther < hubFarther) { kappaRejects++; return; }
                if (farther == hubFarther) tied[remaining++] = tied[i];
            }
            ties = remaining;
        }
        if (REFINEMENT_ROUNDS > 0 && ties > 0) {
            uint64_t colour[MaxOrder];
            refineColours(childOrder, colour);
            const uint64_t hubColour = colour[hub];
            int remaining = 0;
            for (int i = 0; i < ties; ++i) {
                if (colour[tied[i]] < hubColour) { kappaRejects++; return; }
                if (colour[tied[i]] == hubColour) tied[remaining++] = tied[i];
            }
            ties = remaining;
        }

        // -M: exact test
        if constexpr (Maximal) {
            if (childOrder == targetOrder && !isMaximal(childOrder, hubBall, hubSquares)) {
                maximalRejects++;
                return;
            }
        }

        const bool wantGenerators = (childOrder < targetOrder);
        const bool wantCanonical = (ties != 0);
        if (!wantCanonical) shortcutAccepts++;

        if (wantCanonical || wantGenerators)
            runNauty(childOrder, wantCanonical, wantGenerators, hubBall, hubSquares, tied, ties);
        if (wantCanonical) {
            if (orbits[lab[0]] != orbits[hub]) { nautyRejects++; return; }
            nautyAccepts++;
        }

        if (childOrder == targetOrder) writeOut(childOrder);
        else extend(childOrder);
    }

    // =======================================================================
    // The search over feasible pair sets
    // =======================================================================

    void addPair(int hub, int v, int w) {
        partnerAt[hub][v] = w;
        partnerAt[hub][w] = v;
        g.makeAdjacent(hub, v);
        g.makeAdjacent(hub, w);
        g.makeAdjacent(v, w);
        degree[hub] += 2;
        degree[v] += 2;
        degree[w] += 2;
    }

    void removePair(int hub, int v, int w) {
        g.makeNonAdjacent(hub, v);
        g.makeNonAdjacent(hub, w);
        g.makeNonAdjacent(v, w);
        degree[hub] -= 2;
        degree[v] -= 2;
        degree[w] -= 2;
    }

    // Visits every feasible pair set once, in lexicographic order.  'available' are the
    // vertices that can still be chosen.
    void searchPairs(int order, Set available, int k, int limit,
                     Set hubBall, uint64_t hubSquares) {
        // -d#: skip pair sets that are too small or miss a forced vertex
        bool worthTrying = true;
        if constexpr (MinDegree)
            worthTrying = k >= minDegree.leastPairsAt[order]
                       && (minDegree.forcedAt[order] & ~g[order]).isEmpty();
        if (worthTrying) tryChild(order, k, hubBall, hubSquares);
        if (k >= limit) return;

        // -M, last level (plain class): every unmet far triangle needs a chosen vertex of
        // its own.  If there are just enough choices left, all must come from them.
        Set forcedToMeet {};
        if constexpr (Maximal && !C4Free) {
            if (order + 1 == targetOrder) {
                const Set reach = reachOf(g[order], hubBall);
                int unmet = 0;
                for (int i = 0; i < maximal.triangles; ++i) {
                    const Set triangle = maximal.triangle[i];
                    if ((triangle & reach).isNotEmpty()) continue;
                    if ((triangle & available).isEmpty()) return;
                    ++unmet;
                    forcedToMeet |= triangle;
                }
                if (unmet > 2 * (limit - k)) return;
                if (unmet < 2 * (limit - k)) forcedToMeet = Set{};   // still slack
            }
        }

        // -d#: the forced vertices must all be chosen.  If there are just enough choices
        // left, all must come from them.
        Set forcedLeft {};
        if constexpr (MinDegree) {
            forcedLeft = minDegree.forcedAt[order] & ~g[order];
            if ((forcedLeft & ~available).isNotEmpty()) return;
            const int owed = forcedLeft.popcount();
            if (owed > 2 * (limit - k)) return;
            if (owed < 2 * (limit - k)) forcedLeft = Set{};      // still slack
        }

        // A pair set of maximum size must contain every minimum-degree vertex; if that
        // is no longer possible, lower the limit.
        if (limit == maxChooseAt[order]) {
            const Set uncovered = mustCoverAt[order] & ~g[order];
            if ((uncovered & ~available).isNotEmpty()
                || uncovered.popcount() > 2 * (limit - k)) {
                if (k >= --limit) return;
            }
        }

        uint8_t* pairs = pairsAt[order];
        const Graph& c2 = closed2[order];

        // the last pair of a maximum-size set must cover the remaining minimum-degree vertices
        const bool lastPairMustCover = (limit == maxChooseAt[order]) && (limit - k == 1);
        const Set uncovered = lastPairMustCover ? mustCoverAt[order] & ~g[order] : Set{};

        int v;
        while (available.consume_next_bit(v)) {
            Set forbidden;
            if constexpr (C4Free) forbidden = closed3[order][v];
            else                  forbidden = c2[v];

            Set partners = available & ~forbidden;
            if (forcedToMeet.isNotEmpty())
                partners = forcedToMeet.hasBit(v) ? partners & forcedToMeet : Set{};
            if (forcedLeft.isNotEmpty())
                partners = forcedLeft.hasBit(v) ? partners & forcedLeft : Set{};
            if (lastPairMustCover) {
                const Set owed = uncovered & ~Set::bitAt(v);
                if (owed.hasAtLeastWeight<2>()) partners = Set{};   // one w cannot cover two
                else if (owed.isNotEmpty())     partners &= owed;
            }

            for (int w : partners) {
                Set next = available & ~Set::bitAt(w);
                if constexpr (C4Free) next &= ~c2[v] & ~c2[w];
                else                  next &= ~g[v] & ~g[w];

                // the hub's B_2 and sum of squared degrees
                const Set nextHubBall = hubBall | g[v] | g[w] | Set::bitAt(v) | Set::bitAt(w);
                const uint64_t nextHubSquares = hubSquares
                        + uint64_t(degree[v] + 2) * (degree[v] + 2)
                        + uint64_t(degree[w] + 2) * (degree[w] + 2);

                pairs[2 * k] = uint8_t(v);
                pairs[2 * k + 1] = uint8_t(w);
                addPair(order, v, w);
                searchPairs(order, next, k + 1, limit, nextHubBall, nextHubSquares);
                removePair(order, v, w);
            }

            // v is not chosen from here on
            if constexpr (MinDegree)
                if (minDegree.forcedAt[order].hasBit(v)) return;
            if (limit == maxChooseAt[order] && mustCoverAt[order].hasBit(v)) {
                if (k >= --limit) return;
            }
        }
    }

    // =======================================================================
    // One node of the search tree
    // =======================================================================

    void computeClosures(int order) {
        Graph& c2 = closed2[order];
        uint32_t* squares = squaresAt[order];
        for (int v = 0; v < order; ++v) {
            Set reach = g[v] | Set::bitAt(v);
            uint32_t sum = 0;
            for (int w : g[v]) { reach |= g[w]; sum += uint32_t(degree[w]) * degree[w]; }
            c2.adj[v] = reach;
            squares[v] = sum;
        }
        if constexpr (C4Free) {
            Graph& c3 = closed3[order];
            for (int v = 0; v < order; ++v) {
                Set reach = Set::bitAt(v);
                for (int w : g[v]) reach |= c2[w];
                c3.adj[v] = reach;
            }
        }
    }

    // Expands an accepted graph on 'order' vertices; its generators are in
    // generators[order].
    void extend(int order) {
        // res/mod: keep every mod-th node at split.level
        if constexpr (Split)
            if (order == split.level && split.counter++ % split.modulus != split.residue)
                return;

        nodeCount++;
        g.adj[order] = Set{};                 // the hub starts isolated
        degree[order] = 0;

        computeClosures(order);

        // The hub needs minimum degree, so 2k <= delta + 2.  For k = delta/2 + 1 all
        // vertices of degree delta must be chosen.
        int delta = MaxOrder;
        for (int x = 0; x < order; ++x) delta = std::min(delta, degree[x]);

        Set minDegreeSet {};
        for (int x = 0; x < order; ++x) if (degree[x] == delta) minDegreeSet.add(x);

        int maxChoose = order ? delta / 2 + 1 : 0;
        Set mustCover = minDegreeSet;
        if (order && minDegreeSet.popcount() > delta + 2) { --maxChoose; mustCover = Set{}; }
        if (maxChoose > order / 2) { maxChoose = order / 2; mustCover = Set{}; }
        maxChooseAt[order] = maxChoose;
        mustCoverAt[order] = mustCover;

        // -d#: a vertex gains 2 per level when it is chosen, and nothing otherwise
        if constexpr (MinDegree) {
            const int wanted = minDegree.minimum;
            const int levelsLeft = targetOrder - order;

            // (a) a vertex that cannot reach the minimum degree any more
            if (order && delta + 2 * levelsLeft < wanted) return;

            // (b) vertices that must be chosen at this level
            const int critical = wanted - 2 * (levelsLeft - 1);
            Set forced {};
            for (int x = 0; x < order; ++x) if (degree[x] < critical) forced.add(x);
            if (forced.popcount() > 2 * maxChoose) return;

            //     they must be compatible: non-adjacent (C4-free: distance >= 3)
            for (int x : forced) {
                const Set clash = (C4Free ? closed2[order][x] : g[x]) & ~Set::bitAt(x);
                if ((clash & forced).isNotEmpty()) return;
            }

            //     on the last level, each needs a possible partner
            if (levelsLeft == 1)
                for (int x : forced) {
                    const Set reachable = C4Free ? closed3[order][x] : closed2[order][x];
                    if ((Set::bitRange(order) & ~reachable).isEmpty()) return;
                }

            minDegree.forcedAt[order] = forced;
            minDegree.leastPairsAt[order] = critical > 0 ? (critical + 1) / 2 : 0;

            //     the hub itself needs 2k >= critical
            if (minDegree.leastPairsAt[order] > maxChoose) return;

            // (c) Optimistic simulation of the remaining levels: a level with least
            //     degree L chooses at most L+2 vertices; lift the lowest degrees.  The
            //     simulated degrees dominate all reachable ones, so a shortfall here
            //     cannot be repaired.
            int shortfall[2 * MaxOrder + 4] = {};
            int deficient = 0;
            for (int x = 0; x < order; ++x)
                if (degree[x] < wanted) { shortfall[degree[x]]++; ++deficient; }

            int vertices = order;
            for (int i = 0; i < levelsLeft && deficient > 0; ++i) {
                int least = 0;
                while (least < wanted && shortfall[least] == 0) ++least;
                // the cap of this level, as in the computation of maxChoose
                int chosen = (shortfall[least] > least + 2) ? least : least + 2;
                if (i == 0) chosen = 2 * maxChoose;          // known exactly here
                chosen = std::min(chosen, vertices);

                int lifted[2 * MaxOrder + 4] = {};
                int budget = chosen;
                for (int d = least; d < wanted && budget > 0; ++d) {
                    const int move = std::min(budget, shortfall[d]);
                    if (move == 0) continue;
                    shortfall[d] -= move;
                    budget -= move;
                    if (d + 2 < wanted) lifted[d + 2] += move; else deficient -= move;
                }
                for (int d = 0; d < wanted; ++d) shortfall[d] += lifted[d];

                ++vertices;
                if (chosen < wanted) { shortfall[chosen]++; ++deficient; }
            }
            if (deficient > 0) return;
        }

        if constexpr (Maximal)
            if (order + 1 == targetOrder) findFarTriangles(order);

        // degree classes up to delta + 2, the largest possible hub degree
        Set* degreeClass = degreeClassAt[order];
        for (int d = 0; d <= delta + 2; ++d) degreeClass[d] = Set{};
        for (int x = 0; x < order; ++x)
            if (degree[x] <= delta + 2) degreeClass[degree[x]].add(x);
        if (maxChoose > MaxPairsInKey) {          // never reached at any feasible order
            std::cerr << "genll: pair sets of " << maxChoose << " pairs do not fit a key\n";
            exit(1);
        }

        searchPairs(order, Set::bitRange(order), 0, maxChoose, Set{}, 0);
    }

    // the vertices far from v
    Set farFrom(int v, int order) const {
        Set ball = g[v] | Set::bitAt(v);
        for (int r = 1; r < FarRadius; ++r) {
            Set wider = ball;
            for (int w : ball) wider |= g[w];
            ball = wider;
        }
        return Set::bitRange(order) & ~ball;
    }

    // -M: no three vertices of the child are pairwise far
    bool isMaximal(int order, Set hubBall, uint64_t hubSquares) const {
        const int hub = order - 1;
        Set far[MaxOrder], ball;
        for (int v = 0; v < order; ++v) {
            if (v == hub) hubNearInvariant(hub, hubBall, hubSquares, ball);
            else          nearInvariant(hub, v, ball);
            if constexpr (C4Free) { Set wider; ballInvariant(ball, wider); ball = wider; }
            far[v] = Set::bitRange(order) & ~ball;
        }
        for (int v = 0; v < order; ++v)
            for (int w : far[v])
                if (w > v && (far[v] & far[w]).isNotEmpty()) return false;
        return true;
    }

    // -M: a greedy family of vertex-disjoint far triangles of the parent
    void findFarTriangles(int order) {
        Set far[MaxOrder];
        for (int v = 0; v < order; ++v) far[v] = farFrom(v, order);

        maximal.triangles = 0;
        Set unused = Set::bitRange(order);
        for (int v = 0; v < order && maximal.triangles < MaxFarTriangles; ++v) {
            if (!unused.hasBit(v)) continue;
            for (int w : far[v] & unused) {
                const Set third = far[v] & far[w] & unused;
                if (third.isNotEmpty()) {
                    const int x = third.findFirstSetBit();
                    const Set found = Set::bitAt(v) | Set::bitAt(w) | Set::bitAt(x);
                    maximal.corner[maximal.triangles][0] = uint8_t(v);
                    maximal.corner[maximal.triangles][1] = uint8_t(w);
                    maximal.corner[maximal.triangles][2] = uint8_t(x);
                    maximal.triangle[maximal.triangles++] = found;
                    unused &= ~found;
                    break;
                }
            }
        }
    }

    // are x and y no longer far in the child?
    bool closeInChild(int x, int y) const {
        if (g.areAdjacent(x, y)) return true;
        if ((g[x] & g[y]).isNotEmpty()) return true;
        if constexpr (C4Free) {
            Set ball = g[x];
            for (int w : g[x]) ball |= g[w];
            return (ball & g[y]).isNotEmpty();
        }
        return false;
    }

    // A far triangle can only be destroyed if it meets this set: the chosen set, or for
    // C4-free its closed neighbourhood (= hubBall).
    Set reachOf(Set chosen, Set hubBall) const {
        if constexpr (C4Free) return hubBall;
        else                  return chosen;
    }

    void flushOutput() {
        if (outputUsed) { fwrite(outputBuffer, 1, size_t(outputUsed), stdout); outputUsed = 0; }
    }

    void writeOut(int order) {
        ++outCount;
        if constexpr (CountOnly) return;
        if (outputUsed + LongestGraph6 > OutputBufferSize) flushOutput();
        char* p = writeGraph6(g, order, outputBuffer + outputUsed);
        outputUsed = int(p - outputBuffer);
    }
};

// ---------------------------------------------------------------------------

struct Options {
    int order = -1;
    int minDegree = 0;
    int residue = 0, modulus = 1;
    bool c4free = false, c5free = false, countOnly = false, quiet = false, maximal = false;
};

template<bool C4Free, bool Split, bool CountOnly, bool MinDegree, bool Maximal>
static void run(const Options& options) {
    auto generator = std::make_unique<
        Generator<C4Free, Split, CountOnly, MinDegree, Maximal>>();
    generator->targetOrder = options.order;
    generator->c5free = options.c5free;
    if constexpr (MinDegree) generator->minDegree.minimum = options.minDegree;
    if constexpr (Split) {
        generator->split.residue = options.residue;
        generator->split.modulus = options.modulus;
        generator->split.level = std::min(std::max(1, options.order - 3), options.order - 1);
    }

    const auto start = std::chrono::steady_clock::now();
    if (options.order == 0) generator->writeOut(0);
    else                    generator->extend(0);
    generator->flushOutput();
    fflush(stdout);
    const double seconds = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - start).count();

    if (!options.quiet) {
        std::cerr << ">Z " << generator->outCount << " graphs generated in "
                  << seconds << " sec\n";
        if constexpr (GENLL_STATS) {
            std::cerr << ">S nodes=" << generator->nodeCount.get()
                      << " candidates=" << generator->childCount.get()
                      << " nauty=" << generator->nautyCount.get()
                      << " (canonical " << generator->canonCount.get() << ")"
                      << " orbittests=" << generator->orbitTestCount.get()
                      << " largestorbit=" << generator->largestOrbit.get() << "\n";
            const double perCent = 100.0 / double(generator->childCount.get());
            std::cerr << ">K rejected: degree=" << generator->degreeRejects.get() * perCent
                      << "% orbit=" << generator->orbitRejects.get() * perCent
                      << "% kappa=" << generator->kappaRejects.get() * perCent
                      << "% nauty=" << generator->nautyRejects.get() * perCent
                      << "% fartriangle=" << generator->triangleRejects.get() * perCent
                      << "% maximality=" << generator->maximalRejects.get() * perCent
                      << "% | accepted: shortcut=" << generator->shortcutAccepts.get() * perCent
                      << "% nauty=" << generator->nautyAccepts.get() * perCent << "%\n";
        }
    }
}

// The options become template parameters, so unused features cost nothing.
template<bool C4Free, bool Split, bool MinDegree, bool Maximal>
static void runByOutput(const Options& options) {
    if (options.countOnly) run<C4Free, Split, true,  MinDegree, Maximal>(options);
    else                   run<C4Free, Split, false, MinDegree, Maximal>(options);
}
template<bool C4Free, bool Split, bool MinDegree>
static void runByMaximal(const Options& options) {
    if (options.maximal) runByOutput<C4Free, Split, MinDegree, true>(options);
    else                 runByOutput<C4Free, Split, MinDegree, false>(options);
}
template<bool C4Free, bool Split>
static void runByMinDegree(const Options& options) {
    if (options.minDegree > 0) runByMaximal<C4Free, Split, true>(options);
    else                       runByMaximal<C4Free, Split, false>(options);
}
template<bool C4Free>
static void runBySplit(const Options& options) {
    if (options.modulus > 1) runByMinDegree<C4Free, true>(options);
    else                     runByMinDegree<C4Free, false>(options);
}

static void usage() {
    std::cerr <<
        "usage: genll [-u] [-q] [-f] [-5] [-M] [-d#] [res/mod] N   (0 <= N <= " << MaxOrder << ")\n"
        "  -u        count only, write no graphs\n"
        "  -q        no summary on stderr\n"
        "  -f        C4-free locally linear graphs\n"
        "  -5        C5-free locally linear graphs (these are C4-free as well);\n"
        "            not with -M, whose maximality ignores 5-cycles\n"
        "  -M        only maximal ones: no further triangle can be added\n"
        "  -d#       minimum degree at least # (rounded up to an even number)\n"
        "  res/mod   generate only part res of mod; the parts are disjoint and\n"
        "            together complete, and all of them need the same other options\n";
}

int main(int argc, char** argv) {
    Options options;

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        try {
            if      (argument == "-u") options.countOnly = true;
            else if (argument == "-q") options.quiet = true;
            else if (argument == "-f") options.c4free = true;
            else if (argument == "-M") options.maximal = true;
            else if (argument == "-5") options.c5free = true;
            else if (argument.rfind("-d", 0) == 0 && argument.size() > 2)
                options.minDegree = std::stoi(argument.substr(2));
            else if (argument.find('/') != std::string::npos) {
                const size_t slash = argument.find('/');
                options.residue = std::stoi(argument.substr(0, slash));
                options.modulus = std::stoi(argument.substr(slash + 1));
            }
            else if (!argument.empty() && argument[0] != '-')
                options.order = std::stoi(argument);
            else { usage(); return 1; }
        } catch (const std::exception&) { usage(); return 1; }
    }

    if (options.order < 0 || options.order > MaxOrder
        || options.minDegree < 0 || options.minDegree > MaxOrder
        || options.modulus < 1 || options.residue < 0 || options.residue >= options.modulus
        || (options.c5free && options.maximal)) {
        usage();
        return 1;
    }
    if (options.minDegree % 2) ++options.minDegree;   // every degree here is even

    if (!options.quiet) {
        std::cerr << ">A genll";
        if (options.c4free) std::cerr << " -f";
        if (options.c5free) std::cerr << " -5";
        if (options.maximal) std::cerr << " -M";
        if (options.minDegree) std::cerr << " -d" << options.minDegree;
        if (options.modulus > 1) std::cerr << " " << options.residue << "/" << options.modulus;
        std::cerr << " n=" << options.order << "\n";
    }

    if (options.c4free) runBySplit<true>(options);
    else                runBySplit<false>(options);
    return 0;
}
