// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#ifndef TOOLS_H
#define TOOLS_H

#include <cstdint>
#include <assert.h>

//NTL stuff
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_pE.h>
#include <NTL/ZZ_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/tools.h>
using namespace std;
using namespace NTL;

//extending Vec<int64_t>

bool operator< (const Vec<int64_t>& a, const Vec<int64_t>& b);

typedef struct vec_int64_t_less
{
    bool operator()(const Vec<int64_t>& a, const Vec<int64_t>& b) const
    {
        return a < b;
    }
} vi64less;

inline bool operator<(const Vec<int64_t>& a, const Vec<int64_t>& b)
{
    int64_t i, n;
    n = a.length();
    if( n != (int64_t)b.length())
        Error("operator<: dimension mismatch");

    for(i = 0; i < n-1; i++)
    {
        if( a[i] > b[i] )
        {
            return false;
        }
        else if( a[i] < b[i] )
        {
            return true;
        }
    }
    return a[n-1] < b[n-1];
}

template<typename T>
T sum(const Vec<T>& a)
{
    T res = T(0);
    for(int64_t i = 0; i < a.length(); i++)
        res += a[i];
    return res;
}

// res = [0, 1, ..., n - 1] \ B, i.e, the complement of B
// assumes B \subset [0, 1, ..., n - 1]
inline void complement(Vec<int64_t> &res, const int64_t n, const Vec<int64_t> B)
{
    int64_t i, j, position;
    res.SetLength(n - B.length());
    i = 0;
    position = 0;
    for(j = 0; j < B.length(); j++)
    {
        while(i < B[j]  && i < n)
        {
            res[position] = i;
            i++;
            position++;
        }
        i++;
    }
    while(i < n)
    {
        res[position] = i;
        i++;
        position++;
    }
    assert(position == n - B.length());
}

inline void conv(ZZX& x, const zz_pE& a){ conv(x, rep(a)); }
inline void conv(ZZX& x, const ZZ_pE& a){ conv(x, rep(a)); }
inline void conv(zz_pE& x, const ZZX& a){ conv(x, conv<zz_pX>(a)); }
inline void conv(ZZ_pE& x, const ZZX& a){ conv(x, conv<ZZ_pX>(a)); }



#endif
