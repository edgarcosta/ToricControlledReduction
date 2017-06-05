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

template<typename T, typename R, typename S, typename Compare>
inline void conv(map< T, R, Compare> &x, const map< T, S, Compare> &y)
{
    x.clear();
    typename map< T, S, Compare>::const_iterator yit;
    for(yit = y.cbegin(); yit != y.cend(); yit++)
        x[yit->first] = conv<R>(yit->second);
}
}

using namespace NTL;




//returns binomial(n, k)
ZZ binomial(int64_t n, int64_t k);

/*
 * returns factorial(n)/factorial(start-1)
 */
template<typename R>
R  factorial(const int64_t &n, const int64_t &start = 1)
{
    R result=R(1);
    for(int64_t i = start; i <=n ; i++)
        result *= i;
    return result;
}

int64_t valuation_of_factorial(const int64_t n, const int64_t p);

template<typename R> void factorial_padic(R &result, int64_t &val, const int64_t &n, const int64_t &p, const int64_t &start = 1)
{
    result = 1;
    val = 0;
    int64_t i, tmp;
    for(i = start; i <= n; i++)
    {
        if(i%p != 0)
            result *= i;
        else
        {
            tmp = i;
            while(tmp%p == 0)
            {
                tmp /= p;
                val++;
            }
            result *= tmp;
        }
    }
    if( start == 1)
        assert_print(val, ==, valuation_of_factorial(n, p));
}

/*
 * extending Vec<T>
 */
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

template<typename T>
T max(const Vec<T> &v)
{
    T m = v[0];
    for(int64_t i = 1; i < v.length(); i++)
        if(v[i] > m)
            m = v[i];
    return m;
}
template<typename T>
T min(const Vec<T> &v)
{
    T m = v[0];
    for(int64_t i = 1; i < v.length(); i++)
        if(v[i] < m)
            m = v[i];
    return m;
}


template<typename T>
void mul(T &res, const Vec<T> &v, const Vec<T> &w)
{
    assert_print(v.length(), ==, w.length());
    res = T(0);
    T tmp;
    for(int64_t i = 0; i < v.length(); i++)
    {
        mul(tmp, v[i], w[i]);
        add(res, res, tmp);
    }
}
template<typename T>
T operator*(const Vec<T> &v, const Vec<T> &w)
{
    T res;
    mul(res, v, w);
    return res;
}

template<typename T>
void mul(Vec<T> &v, const Mat<T> &A, const Vec<T> &b)
{
    assert_print(b.length(), ==, A.NumCols());
    Vec<T> tmp;
    tmp.SetLength(b.length());
    for(int64_t i = 0; i < b.length(); i++)
        mul(tmp[i], A[i], b);
    swap(tmp, v);
}


template<typename T>
Vec<T> operator*(const Mat<T> &A, const Vec<T> &b)
{
    assert_print(b.length(), ==, A.NumCols());
    Vec<T> v;
    mul(v, A, b); 
    return v;
}

template<typename T>
void add(Vec<T> &res, const Vec<T> &v, const Vec<T> &w)
{
    assert_print(v.length(), ==, w.length());
    res.SetLength(v.length());
    for(int64_t i = 0; i < v.length(); i++)
        add(res[i], v[i], w[i]);
}

template<typename T>
Vec<T> operator+(const Vec<T> &v, const Vec<T> &w)
{
    Vec<T> res;
    add(res, v, w);
    return res;
}

template<typename T>
Vec<T> operator+=(Vec<T> &v, const Vec<T> &w)
{
    add(v, v, w);
    return v;
}

template<typename T>
void sub(Vec<T> &res, const Vec<T> &v, const Vec<T> &w)
{
    assert_print(v.length(), ==, w.length());
    res.SetLength(v.length());
    for(int64_t i = 0; i < v.length(); i++)
        sub(res[i], v[i], w[i]);
}

template<typename T>
Vec<T> operator-(const Vec<T> &v, const Vec<T> &w)
{
    Vec<T> res;
    sub(res, v, w);
    return res;
}

template<typename T>
Vec<T> operator-=(Vec<T> &v, const Vec<T> &w)
{
    sub(v, v, w);
    return v;
}

template<typename T>
void sub(Mat<T> &res, const Mat<T> &A, const Mat<T> &B)
{
    assert_print(A.NumRows(), ==, B.NumRows());
    assert_print(A.NumCols(), ==, B.NumCols());
    res.SetDims(res.NumRows(), res.NumCols());
    for(int64_t i = 0; i < A.NumRows(); i++)
        sub(res[i], A[i], B[i]);
}
template<typename T>
Mat<T> operator-=(Mat<T> &A, const Mat<T> &B)
{
    sub(A, A, B);   
    return A;
}





/*
 * extending << and >>
 */
//istream for a map< Vec<T>, R, Compare>
template<typename T, typename R, typename Compare>
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

//ostream for a map< Vec<T>, R, Compare>
template<typename T, typename R, typename Compare>
NTL_SNS ostream & operator<<(NTL_SNS ostream& s, const  map< Vec<T>, R, Compare>& a)
{
    typename map< Vec<T>, R, Compare>::const_iterator it;
    int64_t i;

    Mat<T>  monomials;
    Vec<R> coefficients;
        
    monomials.SetDims(a.size(), a.cbegin()->first.length());
    coefficients.SetLength(a.size());
    for(i = 0, it = a.cbegin(); it != a.cend(); it++, i++)
    {
        monomials[i] = it->first;
        coefficients[i] = it->second;
    }
    s << monomials<<endl;
    s << coefficients << endl;
    return s;
}

//similar to << but more human readable and compatible with sage 
template<typename T>
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
template<typename T, typename R, typename Compare>
NTL_SNS ostream & operator<<=(NTL_SNS ostream& s, const map< Vec<T>, R, Compare>& a)
{
    s << "{";
    int64_t i, n;
    i = 0;
    n = a.size();
    typename  map< Vec<T>, R, Compare>::const_iterator it;
    for(it = a.cbegin() ; a.cend() != it ; ++it)
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
void complement(Vec<int64_t> &res, const int64_t n, const Vec<int64_t> B);





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
int64_t min_P(const Mat<int64_t> &AP, const Vec<int64_t> &bP, const Vec<int64_t> &v);


//returns the minimal k such that v \in k*int(P)
int64_t min_intP(const Mat<int64_t> &AP, const Vec<int64_t> &bP, const Vec<int64_t> &v);


// reverse dict
template<typename T, typename Compare> void reverse_dict(map<T, int64_t, Compare> &dict, const Vec<T> &v)
{
    for(int64_t i = 0; i < v.length(); i++)
        dict[v[i]] = i;
}


// figure out what ring we are working over
template<typename R> void ring(int64_t &precision, ZZX &fE, ZZ &modulus, const int64_t &p);

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

//lifts it to the respective representative class (ZZ or ZZX) and performs the DivRem operation 
template<typename R> void divrem_lift(R &q, R&r, const R &a, const R&b);

//lifts it to ZZ and performs the DivRem operation
template<typename R> void divrem_lift_ZZ(R &q, R&r, const R &a, const R&b)
{
    ZZ qzz, rzz;
    DivRem(qzz, rzz, conv<ZZ>(a), conv<ZZ>(b));
    q = conv<R>(qzz);
    r = conv<R>(rzz);
}
//lifts it to ZZ and performs the DivRem operation
template<typename R> void divrem_lift_ZZX(R &q, R&r, const R &a, const R&b)
{
    ZZX qzzx, rzzx, azzx;
    DivRem(qzzx, rzzx, conv<ZZX>(a), conv<ZZX>(b));
    q = conv<R>(qzzx);
    r = conv<R>(rzzx);
}

template<> inline void divrem_lift(zz_p &q, zz_p&r, const zz_p &a, const zz_p&b)
{
    divrem_lift_ZZ(q, r, a, b);
}
template<> inline void divrem_lift(ZZ_p &q, ZZ_p&r, const ZZ_p &a, const ZZ_p&b)
{
    divrem_lift_ZZ(q, r, a, b);
}
template<> inline void divrem_lift(zz_pE &q, zz_pE&r, const zz_pE &a, const zz_pE&b)
{
    divrem_lift_ZZX(q, r, a, b);
}
template<> inline void divrem_lift(ZZ_pE &q, ZZ_pE&r, const ZZ_pE &a, const ZZ_pE&b)
{
    divrem_lift_ZZX(q, r, a, b);
}


#endif
