// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#ifndef VEC_INT64_H
#define VEC_INT64_H

#include <cstdint>
#include <assert.h>


#include <NTL/vector.h>
#include <NTL/matrix.h>

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





// res = c * v;
void inline mul(Vec<int64_t> &res, const int64_t c, const Vec<int64_t> &v)
{
    res.SetLength(v.length());
    for(int64_t i = 0; i < v.length(); i++)
        res[i] = c * v[i];
}

// res = v + w
void inline add(Vec<int64_t> &res, const Vec<int64_t> &v, const Vec<int64_t> &w)
{
    assert(v.length() == w.length());
    res.SetLength(v.length());
    for(int64_t i = 0; i < v.length(); i++)
        res[i] = v[i] + w[i];
}
// res = v - w
void inline sub(Vec<int64_t> &res, const Vec<int64_t> &v, const Vec<int64_t> &w)
{
    assert(v.length() == w.length());
    res.SetLength(v.length());
    for(int64_t i = 0; i < v.length(); i++)
        res[i] = v[i] - w[i];
}

// res = v*w
void inline mul(int64_t &res, const Vec<int64_t> &v, const Vec<int64_t> &w)
{
    assert(v.length() == w.length());
    res = 0;
    for(int64_t i = 0; i < v.length(); i++)
        res += v[i]*w[i];
}


// res = A*v
void inline mul(Vec<int64_t> &res, const Mat<int64_t> &A, const Vec<int64_t> &v)
{
    assert(A.NumCols() == v.length());
    res.SetLength(A.NumRows());
    for(int64_t i =0; i < A.NumRows(); i++)
        mul(res[i], A[i], v);
}

/*
 * Operators
 */

inline Vec<int64_t> operator+(const Vec<int64_t>& a, const Vec<int64_t>& b)
{
    Vec<int64_t> x;
    add(x,a,b);
    return x;
}

inline Vec<int64_t> operator-(const Vec<int64_t>& a, const Vec<int64_t>& b)
{
    Vec<int64_t> x;
    sub(x,a,b);
    return x;
}

inline Vec<int64_t> operator*(const int64_t c, const Vec<int64_t>& v)
{
    Vec<int64_t> x;
    mul(x, c , v);
    return x;
}

inline int64_t operator*(const Vec<int64_t> &v, const Vec<int64_t> &w)
{
    int64_t x;
    mul(x, v , w);
    return x;
}

inline Vec<int64_t> operator*(const Mat<int64_t> &A, const Vec<int64_t> &w)
{
    Vec<int64_t> x;
    mul(x, A , w);
    return x;
}




#endif


