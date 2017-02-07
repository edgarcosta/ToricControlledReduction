// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#ifndef TOOLS_H
#define TOOLS_H

#include <cstdint>

//NTL stuff
#include <NTL/matrix.h>
#include <NTL/vector.h>
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

bool operator<(const Vec<int64_t>& a, const Vec<int64_t>& b)
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


#endif
