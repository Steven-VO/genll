// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Steven Van Overberghe
//
// A graph as an array of bitset rows, and its graph6 encoding.
//

#ifndef GENLL_GRAPH_H
#define GENLL_GRAPH_H

#include <array>
#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>

template<typename T>
struct Empty{
    Empty()=default;
    Empty(const Empty &e)=default;
    Empty(T){}

    template<typename S>
    constexpr bool operator==(const S& ){return false;}
};


/*
 * A graph encoded as a list of bitsets.
 * There are three modes:
 *  - max_order=-1 gives graphs on the heap, having arbitrary order
 *  - max_order>0  gives graphs on the stack. Setting exact_order determines whether the max_order is the only permitted order.
 */
template<typename Set, int max_order=-1, bool exact_order= (max_order != -1)>
struct BitGraph{
    typedef int Vertex;

     typedef std::conditional_t<max_order == -1,
        std::vector<Set>,
        std::array<Set, max_order == -1 ? 0: max_order>
    > Arr;
    //typedef ArrayOrVector<Set, order_t> Arr;

    static_assert(max_order >= -1);
    static_assert(max_order != -1 or !exact_order);

    [[no_unique_address]]
    Arr adj;

    [[no_unique_address]]
//    const
    std::conditional_t<exact_order, Empty<int>, int> order_f;

    typedef Set SetType;

    void test_range([[maybe_unused]] Vertex v) const {
        assert(v >= 0 && v < order());
    }

    Set operator[] (Vertex v) const{
        return neighbours(v);
    }

    Set& operator[] (Vertex v) {
        test_range(v);
        return adj[v];
    }

    Set neighbours (Vertex v) const{
        test_range(v);
        return adj[v];
    }

    Set nonNeighbours (Vertex v) const{
        test_range(v);
        return adj[v].complementByXOR(order()) ^ Set::bitAt(v);
    }

    //All neighbours of v that are larger than v
    Set neighbours_up (Vertex v) const {
        test_range(v);
        return adj[v] & Set::fromRange(v);
    }

    //All non-neighbours of v that are larger than v
    Set nonNeighbours_up (Vertex v) const{
        test_range(v);
        return (~adj[v]) & Set::bitRange(v+1, order());
    }

    bool areAdjacent (Vertex v, Vertex w) const {
        test_range(v); test_range(w);
        return adj[v][w];
    }

    //Add an edge between v and w
    void makeAdjacent(Vertex v, Vertex w) {
        test_range(v); test_range(w);
        adj[v].add(w);
        adj[w].add(v);
    }

    void makeNonAdjacent(Vertex v, Vertex w) {
        test_range(v); test_range(w);
        adj[v].remove(w);
        adj[w].remove(v);
    }

    //Create an empty graph of given order
    explicit BitGraph(int order): adj({}), order_f(order) {
        if constexpr (max_order == -1) {
            adj.resize(order);
        }else if constexpr (exact_order) {
            assert(order == max_order);
        }
    }

    BitGraph(): adj({}), order_f(0){}

    constexpr int order() const{
        if constexpr (exact_order){
            return max_order;
        }else{
            return order_f;
        }
    }

    //dynamically change the order of the graph
    void set_order(int order) {
        if constexpr (exact_order) {
            assert(max_order == order);
        }else {
            if constexpr (max_order != -1) {
                assert(0 <= order and order <= max_order);
            }else{
                adj.resize(order);
            }
            order_f = order;
        }
    }

    //number of edges in the graph
    constexpr int size() const{
        int s = 0;
        for (Vertex v = 0; v < order(); ++v) {
            s += adj[v].popcount();
        }
        return s / 2;
    }

    BitGraph complementaryGraph() const{
        BitGraph g (order());
        for (int i = 0; i < order(); ++i) {
            g.adj[i] = nonNeighbours(i);
        }
        return g;
    }

};

// ---------------------------------------------------------------------------
// graph6 output
// ---------------------------------------------------------------------------

// maximum length of the graph6 string (with newline) of a graph on maxOrder vertices
constexpr int longestGraph6(int maxOrder) {
    return 4 + (maxOrder * (maxOrder - 1) / 2 + 5) / 6 + 1;
}

// Writes the graph6 string of the first 'order' vertices of g, plus a newline, at p and
// returns the end.  Needs InverseBitset's bit order: the top x bits of row x are the
// bits graph6 needs, so each row is added with one shift.
template<typename Graph>
inline char* writeGraph6(const Graph& g, int order, char* p) {
    typedef typename Graph::SetType Set;
    constexpr int width = 8 * int(sizeof(Set));
    static_assert(width <= 64, "writeGraph6 needs sets of at most one 64-bit word");
    static_assert((Set::bitAt(0) >> (width - 1)) == 1,
                  "writeGraph6 needs vertex 0 in the top bit of the word, as in InverseBitset");
    constexpr int offset = 63;

    if (order < 63) {
        *p++ = char(order + offset);
    } else {
        *p++ = char(126);
        *p++ = char(((order >> 12) & 63) + offset);
        *p++ = char(((order >> 6) & 63) + offset);
        *p++ = char((order & 63) + offset);
    }

    uint64_t reservoir = 0;      // holds fewer than six bits between rows
    int bits = 0;
    for (int x = 1; x < order; ++x) {
        const uint64_t row = uint64_t(g[x].w) >> (width - x);
        // a long row may need two pieces
        int left = x;
        while (left > 0) {
            const int take = left < 64 - bits ? left : 64 - bits;
            left -= take;
            reservoir = (reservoir << take) | ((row >> left) & ((uint64_t{1} << take) - 1));
            bits += take;
            while (bits >= 6) {
                bits -= 6;
                *p++ = char(((reservoir >> bits) & 63) + offset);
            }
        }
    }
    if (bits) *p++ = char(((reservoir << (6 - bits)) & 63) + offset);
    *p++ = '\n';
    return p;
}

#endif //GENLL_GRAPH_H
