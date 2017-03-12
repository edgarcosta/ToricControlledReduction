// Copyright 2017 Edgar Costa
// See LICENSE file for license details.


#ifndef FINITEDIFF_H
#define FINITEDIFF_H


#include "tools.h"

/*
 * Input:
 * - k, a positive integer > n + 1
 * - G, a vector
 * - M, a vector with (n + 1) matrices, representing  M(Y) = M_0 + M_1 * Y + M_2 Y^2 + ... + M_ * Y^n + M_{n+1} *Y^(n+1)
 * 
 * Output:
 * - H = M(0) M(1) ... M(k-1) G
 *
 *
 * Mainly, works over R, with the goal of having the result correct in S
 */
template<typename R, typename S>
void finitediff(Vec<S> &H, const int64_t &k, const Vec<S> &G, const Vec< Mat<S> > &M)
{
    int64_t l, i, j;
    int64_t n = M.length() - 2;
    int64_t d = G.length();
    assert_print(k, >, n + 1);
    assert_print(G.length(), ==, M[0].NumRows());
    assert_print(G.length(), ==, M[0].NumCols());

    Vec<R> HR;
    HR = conv< Vec<R> >(G);


    // Mfd[l] = M(k - 1 - (self.n + 1) + l) 
    // equivalently
    // Mfd[n + 1 -l] = M(k - 1 - l)
    // Mfd = [ sum( ( self.R(k - 1 - (self.n + 1) + l) ** i) * Mi for i, Mi in  enumerate(M) ) for l in range(self.n + 2) ];
    //
    Vec< Mat<S> > MfdS;
    Vec< Mat<R> > Mfd;
    MfdS.SetLength(n + 2);
    Mfd.SetLength(n + 2);
    for(l = 0; l < n + 2; l++)
    {
        //compute initial table of differences over S
        MfdS[l].SetDims(d, d);
        for(i = 0; i < n + 2; i++)
        {
            S tmp;
            power(tmp,  S(k - 1 - (n + 1) + l), i);
            MfdS[l] += tmp * M[i];
        }
        //then convert them to R
        Mfd[l] = conv< Mat<R> >( MfdS[l] );
    }
    for(l = 0; l < n + 2; l++)
        //u v^{(k - 1 - l)  + 1} --> u v^{(k - 1 - l)}
        // Mfd[n + 1 -l] = M(k - 1 - l)
        HR = Mfd[n + 1 - l] * HR;
    H = conv< Vec<S> >(HR);
    // make Mfd[l] = M[a, a - 1, ..., a - l]
    // where a = k - 1 - (n + 1);
    for(l = 1; l < n + 2; l++)
        for(j = n + 1; j >= l; j--)
            Mfd[j] -= Mfd[j - 1];

    for(l = 0; l < (k - 1 - (n + 1)); l++)
    {
        // Mfd[0] =  M(k - 1 - (n + 1) - l)
        // update Mfd vector
        for(j = n; j >= 0; j--) // deg(M) = n + 1 ==> Mfd[n+1] is constant
            Mfd[j] -= Mfd[j + 1];
        // after
        // Mfd[0] =  M(k - 1 - (n + 1) - l - 1)
        HR = conv< Vec<R> >(H); 
        HR = Mfd[0] * HR;
        H = conv< Vec<S> >(HR);
    }
    //check that we computed M[0] in S
    assert_print( conv< Mat<S> >(Mfd[0]), ==,  M[0]);
}


template<typename S>
void finitediff_lift(Vec<S> &H, const int64_t &k, const Vec<S> &G, const Vec< Mat<S> > &M);

template<>
inline void finitediff_lift(Vec<zz_p> &H, const int64_t &k, const Vec<zz_p> &G, const Vec< Mat<zz_p> > &M)
{
    finitediff<ZZ, zz_p>(H, k, G, M);
}

template<>
inline void finitediff_lift(Vec<ZZ_p> &H, const int64_t &k, const Vec<ZZ_p> &G, const Vec< Mat<ZZ_p> > &M)
{
    finitediff<ZZ, ZZ_p>(H, k, G, M);
}

template<>
inline void finitediff_lift(Vec<zz_pE> &H, const int64_t &k, const Vec<zz_pE> &G, const Vec< Mat<zz_pE> > &M)
{
    finitediff<ZZX, zz_pE>(H, k, G, M);
}

template<>
inline void finitediff_lift(Vec<ZZ_pE> &H, const int64_t &k, const Vec<ZZ_pE> &G, const Vec< Mat<ZZ_pE> > &M)
{
    finitediff<ZZX, ZZ_pE>(H, k, G, M);
}

template<typename R>
void finitediff_plain(Vec<R> &H, const int64_t &k, const Vec<R> &G, const Vec< Mat<R> > &M)
{
    int64_t l, i, j;
    int64_t n = M.length() - 2;
    int64_t d = G.length();
    assert_print(k, >, n + 1);
    assert_print(G.length(), ==, M[0].NumRows());
    assert_print(G.length(), ==, M[0].NumCols());

    H = G;


    // Mfd[l] = M(k - 1 - (self.n + 1) + l) 
    // equivalently
    // Mfd[n + 1 -l] = M(k - 1 - l)
    // Mfd = [ sum( ( self.R(k - 1 - (self.n + 1) + l) ** i) * Mi for i, Mi in  enumerate(M) ) for l in range(self.n + 2) ];
    //
    Vec< Mat<R> > Mfd;
    Mfd.SetLength(n + 2);
    for(l = 0; l < n + 2; l++)
    {
        //compute initial table of differences over S
        Mfd[l].SetDims(d, d);
        for(i = 0; i < n + 2; i++)
        {
            R tmp;
            power(tmp,  R(k - 1 - (n + 1) + l), i);
            Mfd[l] += tmp * M[i];
        }
        //then convert them to R
        Mfd[l] = conv< Mat<R> >( Mfd[l] );
    }
    for(l = 0; l < n + 2; l++)
        //u v^{(k - 1 - l)  + 1} --> u v^{(k - 1 - l)}
        // Mfd[n + 1 -l] = M(k - 1 - l)
        H = Mfd[n + 1 - l] * H;
    
    // make Mfd[l] = M[a, a - 1, ..., a - l]
    // where a = k - 1 - (n + 1);
    for(l = 1; l < n + 2; l++)
        for(j = n + 1; j >= l; j--)
            Mfd[j] -= Mfd[j - 1];

    for(l = 0; l < (k - 1 - (n + 1)); l++)
    {
        // Mfd[0] =  M(k - 1 - (n + 1) - l)
        // update Mfd vector
        for(j = n; j >= 0; j--) // deg(M) = n + 1 ==> Mfd[n+1] is constant
            Mfd[j] -= Mfd[j + 1];
        // after
        // Mfd[0] =  M(k - 1 - (n + 1) - l - 1)
        H = Mfd[0] * H;
    }
    //check that we computed M[0] in S
    assert_print( Mfd[0], ==,  M[0]);
}



#endif
