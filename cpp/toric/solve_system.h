// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#ifndef SOLVE_SYSTEM_H
#define SOLVE_SYSTEM_H


#include "tools.h"

#include <cstdint>
#include <assert.h>

#include <algorithm> // for stl::sort

#include <NTL/ZZX.h>
#include <NTL/vec_long.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/mat_lzz_pE.h>

using namespace std;
using namespace NTL;

//R = ZZ, ZZ_p, zz_p. ZZ_pE, or zz_pE


// if R != ZZ, assumes that R is a field
template<class R>
void pivot_columns(Vec<int64_t> &res, const Mat<R> &T)
{
    //Compute the pivot cols of T
    //doesn't work over R = ZZ
    Mat<R> T_copy = T;
    int64_t rank = gauss(T_copy);
    int64_t ncols = T.NumCols();
    int64_t i, j;
    res.SetLength(rank);
    j = 0;
    for( i = 0; i < rank; i++)
    {
        while(T_copy[i][j] == 0 and j < ncols)
            j++;

        res[i] = j;
        j++;
    }
}

template<> void pivot_columns(Vec<int64_t> &res, const Mat<ZZ> &T);

// X*A = d*I
template<class R>
void inverse(R &d, Mat<R> &X, const Mat<R> &A) 
{
    inv(X, A);
    d = R(1);
}

inline void inverse(ZZ &d, Mat<ZZ> &X, const Mat<ZZ> &A)
{
    inv(d, X, A);
}


/*
 * 
 * R = a domain (Z, Zq or Fq)
 * K = fraction field of R (Q or Qq)
 * T = m * n matriix over R, representing an R-linear map from R^n to R^m
 * 
 * Let T_K = T \otimes_R K, i.e. the corresponding map from K^n to K^m.
 * 
 * This function computes a subset B of {0, 1, ..., m-1}, such that the basis
 * elements e_i of K^m, for i in B, descend to a basis of K^m / im(T_K).
 * 
 * If initB is nonempty then B will be a basis extended from the set initB.
 * 
 * Let J be the matrix of the corresponding inclusion K^B -> K^m, and consider
 * the map T_K + J from K^n \oplus K^B to K^m. By definition of B, this map
 * is surjective.
 * 
 * This function also computes a right inverse of T_K + J, i.e. a map U from
 * K^m to K^n \oplus K^B, such that (T_K + J) U = identity on K^m.
 * 
 * In other words, this function shows how to write every element of R^m
 * as a linear combination of the columns of T, plus possibly some basis
 * elements of R^m not hit by T.
 * 
 * It returns a tuple B, Unum, Uden
 * 
 * Here B is a list containing the indices corresponding to the basis elements,
 * in increasing order.
 * 
 * The matrix U is represented by Unum / Uden, where Unum is a (n + |B|) * m
 * matrix over R, and Uden is a nonzero element of R. The first n rows of U
 * correspond to K^n, and the last |B| rows to K^B.
 */

// R = ZZ, ZZ_p, zz_p. ZZ_pE, or zz_pE
// assumes that p is prime for R = zz_p or zz_pE.
// When p is a prime power one should use solve_system_padic.

template<class R>
void solve_system(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<R> &T, const Vec<int64_t> &initB)
{
    int64_t nrows, ncols;
    int64_t rankS, rankT;
    int64_t i, j;
    Mat<R> S_transpose, S;
    Mat<R> Y; // Y = pivot cols of T, plus J
    Mat<R> Z; // inverse of Y
    Vec<int64_t> pivots_T, pivots_S_transpose;

    nrows = T.NumRows();
    ncols = T.NumCols();

    // deal with a corner case
    if(ncols == 0)
    {
        B.SetLength(nrows);
        clear(Unom);
        Unom.SetDims(nrows,nrows);
        for(i = 0; i < nrows; i++)
        {
            B[i] = i;
            Unom[i][i] = R(1);
        }
        Udenom = R(1);
        return;
    }


    if( initB.length() == 0)
        S = Mat<R>(T);
    else
    {
        // Extends T so that e_i \in Image(T) for i \in initB
        S.SetDims(nrows, ncols + initB.length());
        for(i = 0; i < nrows; i++)
            for(j = 0; j < ncols; j++)
                S[i][j] = T[i][j];
        for(i = 0; i < initB.length(); i++)
            S[ initB[i] ][ ncols + i] = R(1);
    }

    //Compute B
    // B = set of nonpivot rows of S + initB (without swapping rows)
    transpose(S_transpose, S);
    pivot_columns(pivots_S_transpose, S_transpose);
    rankS = pivots_S_transpose.length();

    // B = S.nonpivots() + initB;
    complement(B, nrows, pivots_S_transpose);
    B.append( initB );

    sort(B.begin(), B.end());
    
    for(i = 1; i < B.length(); i++)
        assert( B[i] != B[i-1]);
    
    //Compute the pivot cols of T
    pivot_columns(pivots_T, T);
    rankT = pivots_T.length();
    assert( rankS == rankT + initB.length() ); //otherwise ei \in initB are not LI on K^m/im(T_K)
    assert( rankT + B.length() == nrows );
    

    //Compute part of Y
    // Y = pivot columns of T, plus J
    // where J is the inclusion of B in the R^m 
    Y.SetDims(nrows, nrows);

    //copies pivot columns of T
    for( i = 0; i < nrows; i++)
        for( j = 0; j < rankT; j++)
            //not very cache friendly
            Y[i][j] = T[i][pivots_T[j]];
    // appends the columns corresponding to B
    for( i = 0; i < B.length(); i++)
        Y[B[i]][i + rankT] = R(1);

    // Z = Y.adjoint();
    // Udenom = Y.det()
    inverse(Udenom, Z, Y);
    //Recall
    // m = nrows
    // n = ncols
    // |B| = nrows - rank
    //U from R^m to R^n \oplus R^B, such that (T + J) * Unom = Udenom * identity on R^m
    clear(Unom);
    Unom.SetDims(ncols + B.length(), nrows);
    
    for(i = 0; i < rankT; i++)
        for(j = 0; j < nrows; j++)
            Unom[pivots_T[i]][j] = Z[i][j];

    for(i = rankT; i < nrows; i++)
        for(j = 0; j < nrows; j++ )
            Unom[ncols + i - rankT][j] = Z[i][j];
}

/*
    Same functionality as solve_system(), but first solves over Fq, then lifts
    result to Zq.

    Thus Uden can always be chosen to be a p-adic unit, and in fact this
    function folds it into Unum so that Uden is always 1.
*/
//R = ZZ_p, zz_p. ZZ_pE, or zz_pE, where *_p::modulus = p^precision *_pE::modulus = f
template<class R>
void solve_system_padic(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<R> &T, const Vec<int64_t> &initB, const int64_t p, const int64_t precision, const ZZX f = ZZX(1))
{
    int64_t nrows, ncols;
    int64_t i, j;

    nrows = T.NumRows();
    ncols = T.NumCols();
    if(ncols == 0)
    {
        B.SetLength(nrows);
        Unom.SetDims(nrows,nrows);
        for(i = 0; i < nrows; i++)
        {
            B[i] = i;
            Unom[i][i] = R(1);
        }
        Udenom = R(1);
        return;
    }
    
    Udenom = 1;
    solve_system_local(B, Unom, T, initB, p, f);

    Mat<R> X;
    //Compute X
    X.SetDims(nrows, ncols + B.length() );
    for(i = 0; i < nrows; i++)
        for(j = 0; j < ncols; j++)
            X[i][j] = T[i][j];
    for(i = 0; i < B.length(); i++)
       X[B[i]][i+ncols] = R(1);

    j = 1;
    while( j < precision )
    {
        Unom = 2*Unom - Unom * (X * Unom);
        j *= 2;
    }
}

template<class R>
void solve_system_local(Vec<int64_t> &B,  Mat<R> &Unom, const Mat<R> &T, const Vec<int64_t> &initB, const int64_t p, const ZZX f = ZZX(0))
{
    //f = 0 or f = X
    if( f == 0 or IsX(f))
    {
        //R = ZZ_p or zz_p
        // not an ideal implementation , but I prefer to be careful in case zz_p is being used somewhere else
        Mat<ZZ> T_ZZ = conv< Mat<ZZ> >(T);
        Mat<ZZ> U_ZZ;
        {
            zz_pPush push(p);
            zz_p Udenom;
            Mat<zz_p> U_Fp;
            Mat<zz_p> T_Fp = conv< Mat<zz_p> > (T_ZZ);
            solve_system<zz_p>(B, U_Fp, Udenom, T_Fp, initB);
            assert(Udenom == 1);
            U_ZZ = conv< Mat<ZZ> >(U_Fp);
        }
        Unom = conv< Mat<R> >(U_ZZ);
    }
    else
    {
        //R = ZZ_pE or zz_pE
        Mat<ZZX> T_ZZX;
        int64_t  i, j;
        int64_t ncols, nrows;
        ncols = T.NumCols();
        nrows = T.NumRows();
        T_ZZX = conv< Mat<ZZX> >( T );
        /*
        T_ZZX.SetDims(nrows, ncols);
        for(i = 0; i < nrows; i++)
            for(j = 0; j< ncols; j ++)
                //  T[i][j].rep() \in ZZ_pX or zz_pX
                T_ZZX[i][j] = conv<ZZX>( rep( T[i][j] ));
        */
        Mat<ZZX> U_ZZX;
       {
           zz_pPush push(p);
           {
                zz_pEPush push(conv<zz_pX>(f));
                zz_pE Udenom;
                Mat<zz_pE> U_Fq;
                Mat<zz_pE> T_Fq;
                for(i = 0; i < nrows; i++)
                    for(j = 0; j < ncols; j++)
                        T_Fq[i][j] = conv<zz_pE>( conv<zz_pX>(T_ZZX[i][j]));

                solve_system<zz_pE>(B, U_Fq, Udenom, T_Fq, initB);
                assert(Udenom == 1);
                U_ZZX.SetDims( U_Fq.NumRows(), U_Fq.NumCols());
                for(i = 0; i < U_Fq.NumRows(); i++)
                    for(j = 0; j < U_Fq.NumCols(); j++)
                        U_ZZX[i][j] = conv<ZZX>( rep( U_Fq[i][j] ) );
            }
        }
        Unom = conv< Mat<R> >(U_ZZX);
        /*
        for(i = 0; i < U_ZZX.NumRows(); i++)
            for(j = 0; j < U_ZZX.NumCols(); j++)
                Unom[i][j] = conv< R > ( conv<R::rep_type>( U_ZZX[i][j] ));
        */
    }

}











#endif

