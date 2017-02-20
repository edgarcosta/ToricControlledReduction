// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#ifndef TOOLS_H
#define TOOLS_H

#include <cstdint>
#include <assert.h>
#include <cmath> //ceil
#include <map>
#include <iostream>

#include "vec_int64.h" //extends Vec<int64_t>

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


#ifndef NDEBUG
#define assert_print(left, operator, right) \
{ \
    if( !( (left) operator (right) ) ) \
    { \
        cerr << "ASSERT FAILED: " << #left << " " << #operator << " " << #right << " @ " << __FILE__ << ":" << __LINE__  << endl; \
        cerr << #left << " = " << (left) << "; " << #right << " = " << (right) << endl; \
        abort(); \
    } \
}
#else
#define assert_print(condition, statement) ((void)0)
#endif


#define print(var) { cout << #var << " = " << (var) << endl;}

namespace NTL{
// conversion from ZZX <--> {l}zz_pE
inline void conv(ZZX& x, const zz_pE& a){ conv(x, rep(a)); }
inline void conv(ZZX& x, const ZZ_pE& a){ conv(x, rep(a)); }
inline void conv(zz_pE& x, const ZZX& a){ conv(x, conv<zz_pX>(a)); }
inline void conv(ZZ_pE& x, const ZZX& a){ conv(x, conv<ZZ_pX>(a)); }
}

using namespace NTL;

//extending Vec<T>
template<typename T>
T sum(const Vec<T>& a)
{
    T res = T(0);
    for(int64_t i = 0; i < a.length(); i++)
        res += a[i];
    return res;
}

template<typename T>
T prod(const Vec<T>& a)
{
    T res = T(1);
    for(int64_t i = 0; i < a.length(); i++)
        res *= a[i];
    return res;
}

template<typename T>
Vec<T> operator/(const Vec<T> &v, const T &b)
{
    T i = inv(b);
    return v*i;
}


template<>
inline Vec<ZZ> operator/(const Vec<ZZ> &v, const ZZ &b)
{
    Vec<ZZ> res;
    res.SetLength(v.length());
    for(int64_t i = 0; i < v.length(); i++)
        res[i] = v[i]/b;
    return res;
}


//istream for a map< Vec<T>, R, Compare>
template<class T, class R, class Compare>
NTL_SNS istream & operator>>(NTL_SNS istream& s, map< Vec<T>, R, Compare>& a)
{
    Vec< Vec<T> > monomials;
    Vec<R> coefficients;

    s >> monomials;
    s >> coefficients;
    assert_print(monomials.length(), ==, coefficients.length());
    map< Vec<T>, R, Compare> ibuf;
    for(int64_t i = 0; i < coefficients.length(); i++)
        ibuf[monomials[i]] = coefficients[i];
    a = ibuf;
    return s;
}

//ostream for a map< Vec<uint64_t>, T, vu64less>
template<class T, class R, class Compare>
NTL_SNS ostream & operator<<(NTL_SNS ostream& s, const  map< Vec<T>, R, Compare>& a)
{
    class map< Vec<T>, R, Compare>::const_iterator it;
    int64_t i;

    Vec< Vec<T> > monomials;
    Vec<R> coefficients;
        
    monomials.SetDims(a.size(), a.begin()->first.length());
    coefficients.SetLength(a.size());
    for(i = 0, it = a.begin(); it != a.end(); it++, i++)
    {
        monomials[i] = it->first;
        coefficients[i] = it->second;
    }
    s << monomials<<endl;
    s << coefficients << endl;
    return s;
}

//similar to << but more human readable and compatible with sage 
template<class T>
NTL_SNS ostream & operator<<=(NTL_SNS ostream& s, const Vec<T>& a)
{
    int64_t i, n;
    n = a.length();
    s <<"(";
    for(i = 0; i < n; ++i)
    {
        s << a[i];
        if(i<n-1) s<<", ";
    }
    s << ")";
    return s;
}

//similar to << but more human readable and compatible with sage
template<class T, class R, class Compare>
NTL_SNS ostream & operator<<=(NTL_SNS ostream& s, const map< Vec<T>, R, Compare>& a)
{
    s << "{";
    int64_t i, n;
    i = 0;
    n = a.size();
    typename  map< Vec<T>, R, Compare>::const_iterator it;
    for(it = a.begin() ; a.end() != it ; ++it)
    {
        s <<= it->first;
        s << ": ";
        s << it->second;
        if(i < n - 1)
            s <<",\n ";
        /*else
            break;*/
        i++;
    }
    s << "}";
    return s;
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
    assert_print(position, ==, n - B.length());
}




/*
 * Polyhedron stuff
 */


// Let P be the convex hull of {0} and Pvertices
// st w \in k*P <=>  AP * v + k *bP >= 0 (Half space representation)
//
// returns: points and interior_points vectors of length n st
// points[d] = integral points in d*P 
// interior_points[d] = integral interior points in d*P
// we store them as (n + 1) tuples (d, w) where w \in d*P or equivalently w/d \in P
// d < N
void integral_points(Vec< Vec< Vec< int64_t > > > &points, Vec< Vec< Vec<int64_t> > > &interior_points, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const Vec< Vec<int64_t> > &Pvertices, const int64_t &N);

// returns the minimal k such that v \in k*P
inline int64_t min_P(const Mat<int64_t> &AP, const Vec<int64_t> &bP, const Vec<int64_t> &v)
{
    //v \in k*P iff AP*v + k*bP >= 0;
    int64_t i, j, max, tmp;
    Vec<int64_t> w;
    w.SetLength(bP.length(), 0);
    // v[0] is just a dumb variable to potentially keep track the degree where v came from
    assert_print(AP.NumCols() + 1, ==, v.length());
    assert_print(AP.NumRows(), ==,  bP.length());
    // w = Ap * v[1:]
    for(i = 0; i < AP.NumRows(); i++)
        for(j = 0; j < AP.NumCols(); j++)
            w[i] += AP[i][j] * v[j + 1];
    max = -1;
    for(i = 0; i < bP.length(); i++)
    {
        if( w[i] < 0 and bP[i] == 0)
            return INT64_MAX;//a practical +Infinity

        if( bP[i] != 0)
        {
            tmp = ceil(static_cast<double>( -w[i] ) / bP[i] );
            if( tmp > max )
                max = tmp;
        }
    }
    return max;
}

//returns the minimal k such that v \in k*int(P)
inline int64_t min_intP(const Mat<int64_t> &AP, const Vec<int64_t> &bP, const Vec<int64_t> &v)
{
    //v \in k*int(P) iff AP*v + k*bP > 0;
    int64_t i, j, max, tmp;
    Vec<int64_t> w;
    w.SetLength(bP.length(), 0);
    // v[0] is just a dumb variable to potentially keep track the degree where v came from
    assert_print(AP.NumCols() + 1, ==, v.length());
    assert_print(AP.NumRows(), ==,  bP.length());
    // w = Ap * v[1:]
    for(i = 0; i < AP.NumRows(); i++)
        for(j = 0; j < AP.NumCols(); j++)
            w[i] += AP[i][j] * v[j + 1];
    
    max = -1;
    for(i = 0; i < bP.length(); i++)
    {
        if( (w[i] < 1) and (bP[i] == 0) )
            return INT64_MAX;//a practical +Infinity

        if( bP[i] != 0)
        {
            tmp = ceil(static_cast<double>( -(w[i] - 1) ) / bP[i] );
            if( tmp > max )
                max = tmp;
        }
    }
    return max;
}

// reverse dict
template<class T, class Compare> void reverse_dict(map<T, int64_t, Compare> &dict, const Vec<T> &v)
{
    for(int64_t i = 0; i < v.length(); i++)
        dict[v[i]] = i;
}


// figure out what ring we are working over
template<class R> void ring(int64_t &precision, ZZX &fE, ZZ &modulus, const int64_t &p);

template <> inline void ring<ZZ>(int64_t &precision, ZZX &fE, ZZ &modulus, const int64_t &p)
{
    assert_print( p, ==, 0);
    precision = 0;
    fE = 0;
    modulus = 0;
}


template <> inline void ring<zz_p>(int64_t &precision, ZZX &fE, ZZ &modulus, const int64_t &p)
{
    assert_print( p, !=, 0);
    fE = 0;
    modulus = zz_p::modulus();
    precision = llround( log(modulus)/log(p) );

}


template <> inline void ring<ZZ_p>(int64_t &precision, ZZX &fE, ZZ &modulus, const int64_t &p)
{
    assert_print( p, !=, 0);
    fE = 0;
    modulus = ZZ_p::modulus();
    precision = llround( log(modulus)/log(p) );
}



template <> inline void ring<zz_pE>(int64_t &precision, ZZX &fE, ZZ &modulus, const int64_t &p)
{
    assert_print( p, !=, 0);
    fE = conv<ZZX>(zz_pE::modulus());
    modulus = zz_p::modulus();
    precision = llround( log(modulus)/log(p) );

}

template <> inline void ring<ZZ_pE>(int64_t &precision, ZZX &fE, ZZ &modulus, const int64_t &p)
{
    assert_print( p, !=, 0);
    fE = conv<ZZX>(ZZ_pE::modulus());
    modulus = ZZ_p::modulus();
    precision = llround( log(modulus)/log(p) );

}





#endif
