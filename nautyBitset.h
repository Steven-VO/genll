// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) -2026 Steven Van Overberghe
//
// A bitset of one machine word with nauty's bit order: element i is bit width-1-i,
// so an array of these is a nauty graph as it stands.
//

#ifndef INVERSEBITSET_H
#define INVERSEBITSET_H


#include <cstdint>
#include <array>
#include <type_traits>
#if defined(__LZCNT__)
#include <immintrin.h>
#endif
#include <vector>
#include <ostream>
#include <cassert>

template<int maxSize>
struct InverseBitsetIterator;

template<int maxSize=64>
struct InverseBitset{
    typedef std::conditional_t<maxSize == 32, uint32_t,
            std::conditional_t<maxSize == 64, uint64_t, unsigned __int128>> bitword;

    [[no_unique_address]]
    bitword w;

    static_assert(maxSize == 32 or maxSize == 64 or maxSize == 128);

    constexpr static bitword bitAt(int i){
        return bitword{1} << (maxSize - 1 - i);
    }

    //Construct a bitset of which the first 'length' bits are set
    constexpr static InverseBitset bitRange(int length){
        assert(length >= 0 and length <= maxSize);
        if (length == 0){return InverseBitset();}
        return (~InverseBitset{}) << (maxSize - length);
    }

    //Construct a bitset of which a range of bits is set
    constexpr static InverseBitset bitRange(int startIncl, int stopExcl){
        assert(stopExcl >= startIncl);
        return bitRange(stopExcl - startIncl) >> startIncl;
    }

    //Construct a bitset for which all bits starting from 'startIncl' and higher are set
    constexpr static InverseBitset fromRange(int startIncl){
        return (~InverseBitset{}) >> startIncl;
    }

    constexpr bool hasBit(int i) const{
        return (w & bitAt(i)) != bitword{0};
    }

    constexpr bool operator[](int i) const{
        return hasBit(i);
    }

    void addBit(int i){
        w |= bitAt(i);
    }

    void add(int i){
        addBit(i);
    }

    void removeBit(int i){
        w &= ~(bitAt(i));
    }

    void remove(int i){
        removeBit(i);
    }

    void toggleBit(int i){
        w ^= bitAt(i);
    }

    //Return the bits in range [0..n[ not present in this bitset
    InverseBitset complement(int n) const{
        return bitRange(n) & (~w);
    }

    //If you know this set is a subset of [0..n[
    InverseBitset complementByXOR(int n) const{
        return bitRange(n) ^ (*this);
    }

    constexpr InverseBitset without(int i) const{
        return InverseBitset(w & ~bitAt(i));
    }

    //--------weights------------------

    constexpr bool isEmpty() const{
        return w == bitword{0} ;
    }

    constexpr bool isNotEmpty() const{
        return !isEmpty();
    }

    constexpr int popcount() const{
        if constexpr(sizeof(bitword) == 16){
        return __builtin_popcountll(w) + __builtin_popcountll(w >> 64);;
        }else if constexpr(sizeof(bitword) == 4) {
            return __builtin_popcount(w);
        } else {
            return __builtin_popcountll(w);
        }
    }

    constexpr int size() const{
        return popcount();
    }

    /*
     * Check if the number of set bits in w is equal to a predefined number
     */
    constexpr bool hasWeight(int weight) const{
        return popcount() == weight;
    }

    /*
     * Check if the number of set bits in w is equal to a predefined constexpr number
     */
    template<int weight>
    constexpr bool hasWeight() const{
        static_assert(weight >= 0);

        if constexpr (weight == 0){
            return isEmpty();
        }else if constexpr (weight == 1){
            return (!(w & (w-1))) && w;
    //        return w && !(w & (w-1));
        }else{
            return hasWeight(weight);
        }
    }

    /*
     * Check if the number of set bits in w is at least a predefined number
     */
    template<int weight>
    constexpr bool hasAtLeastWeight() const{
        static_assert(weight >= 0);

        if constexpr (weight <= 0){
            return true;
        }else if constexpr (weight == 1) {
            return isNotEmpty();
        }else if constexpr (weight == 2) {
            return withoutLowest().isNotEmpty();
        }else{
            return popcount() >= weight;
        }
    }

    //Find the index of the first set bit in this bitset.  Undefined if the bitset is empty.
    constexpr int findFirstSetBit() const{
        assert(isNotEmpty());
        if constexpr (maxSize == 32) {
            return __builtin_clz(w);
        }else if constexpr (maxSize == 64) {
            return __builtin_clzll(w);
        }else{
            uint64_t high = uint64_t(w >> 64);
            return high ? __builtin_clzll(high) : 64 + __builtin_clzll(uint64_t(w));
        }
    }

    //Find the index of the first set bit, or maxSize if the bitset is empty (as lzcnt does)
    int findFirstSetBitOrEnd() const{
#if defined(__LZCNT__)
        if constexpr (maxSize == 32) {
            return _lzcnt_u32(w);
        }else if constexpr (maxSize == 64) {
            return _lzcnt_u64(w);
        }
#endif
        return isEmpty() ? maxSize : findFirstSetBit();
    }

    //Return this bitset but with the lowest-order bit cleared
    constexpr InverseBitset withoutLowest() const{
        if (isEmpty()){return bitword{0};}
        return w ^ bitAt(findFirstSetBit());
    }

    //Return the lowest-order bit (as a bitset) present in this bitset
    constexpr InverseBitset lowestOrderBit() const{
        return w ^ withoutLowest();
    }

    //Find first set bit, toggle it, and set argument to the index of this set bit.
    //Return whether this bitset was non-empty
    // if 'safe', prevents undefined (but ignored) behaviour at the cost of an extra branch
    template<bool safe = true>
    bool consume_next_bit(int& index){
        bool cont = isNotEmpty();
        if constexpr (safe){
            if(!cont){return false;}
        }
        index = findFirstSetBit();
        toggleBit(index);
        return cont;
    }

    constexpr InverseBitsetIterator<maxSize> begin() const;
    constexpr static InverseBitsetIterator<maxSize> end() ;

    //-----------Constructors-------------

    //Create a bitword from an array of the indices of all present bits
    template<typename T, unsigned long size, int i=size>
    constexpr static bitword bitWordFromArray(const std::array<T, size>& indices){
        static_assert(i >= 0);

        if constexpr (i == 0){
            return bitword{0};
        }else{
            return bitAt(indices[i-1]) | bitWordFromArray<T,size,i-1>(indices);
        }
    }


    constexpr InverseBitset(const bitword& word): w(word) {}
    constexpr InverseBitset(): w(bitword{0}){}

    explicit InverseBitset(const std::vector<int>& indices){
        w = (bitword) 0;
        for(int i: indices){w |= bitAt(i);}
    }

    template<typename T, unsigned long size>
    constexpr explicit InverseBitset(const std::array<T, size>& indices): w(bitWordFromArray(indices)) {}



    constexpr InverseBitset operator& (const InverseBitset& other) const{
        return InverseBitset(w & other.w);
    }

    constexpr void operator &= (const InverseBitset& other){
        w &= other.w;
    }

    constexpr InverseBitset operator| (const InverseBitset& other) const{
        return InverseBitset(w | other.w);
    }

    constexpr void operator |= (const InverseBitset& other){
        w |= other.w;
    }

    constexpr InverseBitset operator ^ (const InverseBitset& other) const{
        return InverseBitset(w ^ other.w);
    }

    constexpr void operator ^= (const InverseBitset& other){
        w ^= other.w;
    }

    constexpr InverseBitset operator~() const{
        return InverseBitset(~w);
    };

    constexpr bool operator==(const InverseBitset& other) const{
        return w == other.w;
    }

    constexpr bool operator!=(const InverseBitset& other) const{
        return w != other.w;
    }

    constexpr InverseBitset operator<< (int x) const{
        return w << x;
    }

    constexpr InverseBitset operator>> (int x) const{
        return w >> x;
    }

//    constexpr operator unsigned long () const{
//        return w;
//    }

//    constexpr operator bool() const{
//        return !isEmpty();
//    }


    //Make a bitset printable
    friend std::ostream& operator<<(std::ostream& os, const InverseBitset& set){
        for (int i = 0; i < maxSize; ++i) {
            if(set[i]){
                os << "1";
            }else{
                os << "0";
            }
        }
        return os;
    }
};

//Iterates over the elements; the end is index maxSize
template<int maxSize=64>
struct InverseBitsetIterator{
    InverseBitset<maxSize> current;
    int index;

    explicit InverseBitsetIterator(const InverseBitset<maxSize>& set):
        current(set), index(set.findFirstSetBitOrEnd()) {}

    constexpr InverseBitsetIterator(): current(), index(maxSize) {}

    constexpr bool operator != (const InverseBitsetIterator& other) const{
        return index != other.index;
    }

    void operator++(){
        current.toggleBit(index);
        index = current.findFirstSetBitOrEnd();
    }

    constexpr int operator*() const{
        return index;
    }
};

template<int maxSize>
constexpr InverseBitsetIterator<maxSize> InverseBitset<maxSize>::begin() const{
    return InverseBitsetIterator<maxSize>{*this};
}

template<int maxSize>
constexpr InverseBitsetIterator<maxSize> InverseBitset<maxSize>::end() {
    return InverseBitsetIterator<maxSize>{};
}

#endif
