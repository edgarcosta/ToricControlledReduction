// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H


#include "tools.h"

#include <cstdint>
#include <assert.h>

#include <algorithm> // for std::sort and std::binary_search

#include <NTL/LLL.h>
#include <NTL/ZZX.h>
#include <NTL/vec_long.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/mat_lzz_pE.h>

using namespace std;
using namespace NTL;

//extending Mat<ZZ>
//
//computes the characteristic polynomial of M using Newton identities
// = det(I*X - M) = X^n - trace(M)*X^(n - 1) + ...
Vec<ZZ> charpoly(const Mat<ZZ> &M);
// computes the trace
ZZ trace(const Mat<ZZ> &M);


// works for R = ZZ_p and zz_p
template<typename R>
int64_t find_pivot_modp(const Mat<R> &A, const int64_t start_row, const int64_t end_row, const int64_t col, const int64_t p)
{
    for(int64_t r = start_row; r < end_row; ++r)
      if( not is_zero_modp(A[r][col], p) )
        return r;

    return -1;
}


// R is a p-adic ring
// Performs unitary row operations so as to bring A into row echelon
// form (mod p).
// returns the rank
template<typename R>
int64_t gauss_padic(Mat<R> &A, const int64_t p)
{

  int64_t m = A.NumRows();
  int64_t n = A.NumCols();
  int64_t rank = 0;
  int64_t pivot_row = 0;
  int64_t pivot_col = 0;

  while( pivot_row < m and pivot_col < n) {
    int64_t r = find_pivot_modp(A, pivot_row, m, pivot_col, p);

    if(r == -1) {
      pivot_col++;
      continue;
    }

    if (r != pivot_row)
      swap(A[r], A[pivot_row]);

    rank++;

    // rescale the pivot row
    R pivot_inverse = 1/A[pivot_row][pivot_col];
    
    A[pivot_row][pivot_col] = 1;
    for(int64_t j = pivot_col + 1; j < n; ++j)
      A[pivot_row][j] *= pivot_inverse;

    for(int64_t i = pivot_row + 1; i < m; ++i)
    {
      //A[i] -= A[pivot_row] * A[i][pivot_col]
      for(int64_t j = pivot_col + 1; j < n; ++j)
          A[i][j] -= A[pivot_row][j]*A[i][pivot_col];

      A[i][pivot_col] = 0;
    }

    pivot_row++;
    pivot_col++;
  }
  return rank;
}


// R is a p-adic ring
// Performs unitary row operations on A to bring A into row echelon
// form (mod p).
// L is the matrix such that A_new = L * A_old
// the function returns the rank of A
template<typename R>
int64_t LU_padic(Mat<R> &L, Mat<R> &A, const int64_t p, const Vec<int64_t> &skip)
{

  Mat<R> A_copy = A;
  int64_t m = A.NumRows();
  int64_t n = A.NumCols();
  int64_t rank = 0;
  int64_t pivot_row = 0;
  int64_t pivot_col = 0;

  clear(L);
  L.SetDims(m, m);
  for(int64_t i = 0; i < m; ++i)
    L[i][i] = 1;

  Mat<R> LA;
  Mat<R> Aold, Lold;
  while( pivot_row < m and pivot_col < n) {
    if(binary_search(skip.begin(), skip.end(), pivot_col)) {
        pivot_col++;
        continue;
    }
    int64_t r = find_pivot_modp(A, pivot_row, m, pivot_col, p);

    if(r == -1) {
      pivot_col++;
      continue;
    }


    if (r != pivot_row) {
      swap(A[r], A[pivot_row]);
      swap(L[r], L[pivot_row]);
    }

    /*
    Aold = A; Lold = L;
    LA = L * A_copy;
    if( LA != A)
    {
      cout<<endl<<endl<<endl;
      cout<<"Before inverse"<<endl;
      print(pivot_row);
          print(pivot_col);
          print_raw(A_copy);
          print_raw(LA);
          print_raw(L);
          print_raw(A);
          print_raw(Lold);
          print_raw(Aold);
      abort();
    }
    */


    rank++;

    // rescale the pivot row
    R pivot_inverse = 1/A[pivot_row][pivot_col];

    //A[pivot_row][pivot_col] = 1;
    //A[pivot_row][pivot_col] *= pivot_inverse;
    for(int64_t j = 0; j < n; ++j)
      A[pivot_row][j] *= pivot_inverse;
    /*
    if( skip.length() > 0)
      for(int64_t j = 0; j < skip.length() and skip[j] < pivot_col; ++j)
         A[pivot_row][skip[j]] *= pivot_inverse;
         */

    for(int64_t j = 0; j < m; ++j)
      L[pivot_row][j] *= pivot_inverse;

    /*
    LA = L * A_copy;
    if( LA != A)
    {
      cout<<endl<<endl<<endl;
      cout<<"After inverse"<<endl;
      print(pivot_row);
          print(pivot_col);
          print_raw(A_copy);
          print_raw(LA);
          print_raw(L);
          print_raw(A);
          print_raw(Lold);
          print_raw(Aold);
      abort();
    }*/

    for(int64_t i = 0; i < m; ++i) {
      //FIXME
      //Aold = A; Lold = L;
      if( i != pivot_row ) {
        //A[i] -= A[pivot_row] * A[i][pivot_col]
        R c = A[i][pivot_col];
        for(int64_t j = 0; j < n; ++j)
            A[i][j] -= A[pivot_row][j]*c;// - A[i][j];

        /*
        if( skip.length() > 0)
          for(int64_t j = 0; j < skip.length() and skip[j] < pivot_col; ++j)
            A[i][skip[j]] = A[pivot_row][skip[j]]*A[i][pivot_col] - A[i][skip[j]];
            */


        //L[i] -= L[pivot_row] * A[i][pivot_col]
        for(int64_t j = 0; j < m; ++j)
          L[i][j] -= L[pivot_row][j]*c;// - L[i][j];

       /* 
        LA = L * A_copy;
        if( LA != A)
        {
          cout<<endl<<endl<<endl;
          cout<<"After row = "<<i<<endl;
          print(pivot_row);
          print(pivot_col);
          print_raw(A_copy);
          print_raw(LA);
          print_raw(L);
          print_raw(A);
          print_raw(Lold);
          print_raw(Aold);
          abort();
        }
        */

      }
      /*
      //double check here!
      LA = L * A_copy;
      if( LA != A)
      {
        cout<<endl<<endl<<endl;
        cout<<"at the end"<<endl;
        print_raw(LA);
        print_raw(L);
        print_raw(A);
        abort();
      }
      */
    }
    pivot_row++;
    pivot_col++;
  }

  assert_print(L * A_copy, ==, A);
  return rank;
}

//returns the A.NumRows() or 0, if not full rank
template<typename R>
int64_t inverse_padic(Mat<R> &Ainv, Mat<R> A, const int64_t p)
{
  int64_t m = A.NumRows();
  int64_t n = A.NumCols();
  assert_print(m, ==, n);
  int64_t rank = 0;
  int64_t pivot_row = 0;
  int64_t pivot_col = 0;

  clear(Ainv);
  Ainv.SetDims(m, n);
  for(int64_t i = 0; i < m; ++i)
    Ainv[i][i] = 1;

  while( pivot_row < m and pivot_col < n) {
    int64_t r = find_pivot_modp(A, pivot_row, m, pivot_col, p);

    if(r == -1)
      return 0;

    if (r != pivot_row) {
      swap(A[r], A[pivot_row]);
      swap(Ainv[r], Ainv[pivot_row]);
    }

    rank++;


    // rescale the pivot row
    R pivot_inverse = 1/A[pivot_row][pivot_col];
    A[pivot_row][pivot_col] = 1;
    for(int64_t j = pivot_col + 1; j < n; ++j)
      A[pivot_row][j] *= pivot_inverse;
    for(int64_t j = 0; j < n; ++j)
      Ainv[pivot_row][j] *= pivot_inverse;

    // Forward and back substitution
    for(int64_t i = 0; i < m; ++i) {
      if(i != pivot_row) {
        //A[i] -= A[pivot_row] * A[i][pivot_col]
        for(int64_t j = pivot_col + 1; j < n; ++j)
          A[i][j] -= A[pivot_row][j]*A[i][pivot_col];
        //Ainv[i] -= Ainv[pivot_row] * A[i][pivot_col]
        for(int64_t j = 0; j < n; ++j)
          Ainv[i][j] -= Ainv[pivot_row][j]*A[i][pivot_col];
        A[i][pivot_col] = 0;
      }
    }
    pivot_row++;
    pivot_col++;
  }
  assert_print(rank, ==, n);
  return rank;
}

//R = ZZ, ZZ_p, zz_p. ZZ_pE, or zz_pE
// if R != ZZ, assumes that R is a field
template<typename R>
void pivot_columns(Vec<int64_t> &res, const Mat<R> &T) {
    //Compute the pivot cols of T
    //doesn't work over R = ZZ, see specialization
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

// Specialization for R = ZZ
template<> void pivot_columns(Vec<int64_t> &res, const Mat<ZZ> &T);

//R = ZZ, ZZ_p, zz_p. ZZ_pE, or zz_pE
// Assumes that T is row reduced
template<typename R>
void pivot_columns_rr(Vec<int64_t> &res, const Mat<R> &T, const Vec<int64_t> &skip) {
    int64_t ncols = T.NumCols();
    int64_t nrows = T.NumRows();
    int64_t j = 0;
    int64_t k = 0;
    res.SetLength(min(ncols, nrows));
    for(int64_t i = 0; i < nrows and j < ncols; i++) {
        while((binary_search(skip.begin(), skip.end(), j) or T[i][j] == 0) and j < ncols)
            j++;

        if( j < ncols ) {
          res[k] = j;
          k++;
        }
        j++;
    }
    res.SetLength(k);
}


// R is a p-adic ring, and returns the pivot columns modp
template<typename R>
void pivot_columns_modp(Vec<int64_t> &res, const Mat<R> &T, const int64_t p) {
    //Compute the pivot cols of T
    Mat<R> T_copy = T;
    int64_t rank = gauss_padic(T_copy, p);
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

// X*A = d*I
template<typename R>
void inverse(R &d, Mat<R> &X, const Mat<R> &A) {
    inv(X, A);
    d = R(1);
}
template<>
inline void inverse(ZZ &d, Mat<ZZ> &X, const Mat<ZZ> &A)
{
    inv(d, X, A);
}

// K = T.left_kernel();
inline void kernel(Mat<ZZ> &K,  Mat<ZZ> &T) {
    ZZ det2;
    int64_t r = image(det2, T, K);
    K.SetDims(K.NumRows() - r, K.NumCols());
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

template<typename R>
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

    // Compute B
    // B = set of nonpivot rows of S + initB (without swapping rows)
    transpose(S_transpose, S);
    pivot_columns(pivots_S_transpose, S_transpose);
    rankS = pivots_S_transpose.length();

    // B = S.nonpivots() + initB;
    complement(B, nrows, pivots_S_transpose);
    B.append( initB );

    sort(B.begin(), B.end());

    for(i = 1; i < B.length(); i++)
        assert_print( B[i], !=, B[i-1]);

    //Compute the pivot cols of T
    pivot_columns(pivots_T, T);
    rankT = pivots_T.length();
    assert_print( rankS, ==, rankT + initB.length() ); //otherwise ei \in initB are not LI on K^m/im(T_K)
    assert_print( rankT + B.length(), ==, nrows );


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
 Same functionally as solve_system(), but over a p-adic ring
 */

//FIXME
template<typename R>
void solve_system_padic2(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<R> &T,  Vec<int64_t> initB, const int64_t p, const int64_t precision, const ZZX f = ZZX(0))
{
    assert_print(precision, >, 0);
    print(precision);
    //f = 0 or f = X
    assert_print( f == 0, or, IsX(f));
    int64_t nrows = T.NumRows();
    int64_t ncols = T.NumCols();
    Mat<R> Y; // Y = pivot cols of T, plus J
    Mat<R> Z; // inverse of Y
    Vec<int64_t> pivots_T_transpose;


    // deal with a corner case
    if(ncols == 0)
    {
        B.SetLength(nrows);
        clear(Unom);
        Unom.SetDims(nrows,nrows);
        for(int64_t i = 0; i < nrows; i++)
        {
            B[i] = i;
            Unom[i][i] = R(1);
        }
        Udenom = R(1);
        return;
    }

    sort(initB.begin(), initB.end());

    //Compute B
    // B = set of nonpivot rows of S + initB (without swapping rows)
    Mat<R> T_transpose;
    transpose(T_transpose, T);
    Mat<R> PL;
    // Row reduce S_transpose
    int64_t rankT = LU_padic(PL, T_transpose, p, initB);
    // pivots_T_transpose \subset complement(initB)
    // as we forced the columns in initB to not be pivot
    pivot_columns_rr(pivots_T_transpose, T_transpose, initB);

    assert_print(rankT, ==, pivots_T_transpose.length());

    // B = T.nonpivots() + initB;
    complement(B, nrows, pivots_T_transpose);
    sort(B.begin(), B.end());


    clear(Unom);
    Unom.SetDims(ncols + B.length(), nrows);
    Udenom = 1;


    for(int64_t i = 0; i < B.length(); ++i)
      Unom[ncols + i][B[i]] = 1;


    for(int64_t i = 0; i < pivots_T_transpose.length(); ++i) {
      if( ! binary_search(initB.begin(), initB.end(), pivots_T_transpose[i]) ) {

        for(int64_t j = 0; j < ncols; ++j)
          Unom[j][pivots_T_transpose[i]] = PL[i][j];

        for(int64_t j = 0; j < B.length(); ++j)
          Unom[ncols + j][pivots_T_transpose[i]] = -T_transpose[i][B[j]];
      }
    }


    Mat<R> newS;

    newS.SetDims(nrows, ncols + B.length());
    for(int64_t i = 0; i < nrows; i++)
        for(int64_t j = 0; j < ncols; j++)
            newS[i][j] = T[i][j];
    for(int64_t i = 0; i < B.length(); i++)
        newS[ B[i] ][ ncols + i] = R(1);
    print(B);
    print(initB);
    Mat<R> eye;
    eye.SetDims(nrows, nrows);
    for(int64_t i = 0; i < nrows; i++)
      eye[i][i] = 1;

    assert_print(newS * Unom, ==, eye);

}

/*
    Same functionality as solve_system(), but first solves over Fq, then lifts
    result to Zq.

    Thus Uden can always be chosen to be a p-adic unit, and in fact this
    function folds it into Unum so that Uden is always 1.
*/
//R = ZZ_p, zz_p. ZZ_pE, or zz_pE, where *_p::modulus = p^precision *_pE::modulus = f

template<typename R>
void solve_system_local(Vec<int64_t> &B,  Mat<R> &Unom, const Mat<R> &T, const Vec<int64_t> &initB, const int64_t p, const ZZX f = ZZX(0));

template<typename R>
void solve_system_padic(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<R> &T, const Vec<int64_t> &initB, const int64_t p, const int64_t precision, const ZZX f = ZZX(0))
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


//R = ZZ_p or zz_p
template<typename R>
void solve_system_local_lzz_p(Vec<int64_t> &B,  Mat<R> &Unom, const Mat<R> &T, const Vec<int64_t> &initB, const int64_t p, const ZZX f = ZZX(0))
{
    //f = 0 or f = X
    assert_print( f == 0, or, IsX(f));
    // not an ideal implementation , but I prefer to be careful in case zz_p is being used somewhere else
    Mat<ZZ> T_ZZ = conv< Mat<ZZ> >(T);
    Mat<ZZ> U_ZZ;
    {
        zz_pPush push(p);
        zz_p Udenom;
        Mat<zz_p> U_Fp;
        Mat<zz_p> T_Fp = conv< Mat<zz_p> > (T_ZZ);
        solve_system<zz_p>(B, U_Fp, Udenom, T_Fp, initB);
        assert_print(Udenom, ==, 1);
        U_ZZ = conv< Mat<ZZ> >(U_Fp);
    }
    Unom = conv< Mat<R> >(U_ZZ);
}

//R = ZZ_pE or zz_pE
template<typename R>
void solve_system_local_lzz_pE(Vec<int64_t> &B,  Mat<R> &Unom, const Mat<R> &T, const Vec<int64_t> &initB, const int64_t p, const ZZX f = ZZX(0))
{
    assert_print( f, !=, 0);
    assert_print( IsX(f), ==, false);

    Mat<ZZX> T_ZZX;
    T_ZZX = conv< Mat<ZZX> >( T );
    Mat<ZZX> U_ZZX;
   {
       zz_pPush push1(p);
       {
            zz_pEPush push2(conv<zz_pX>(f));
            zz_pE Udenom;
            Mat<zz_pE> U_Fq;
            Mat<zz_pE> T_Fq;
            T_Fq = conv< Mat<zz_pE> >(T_ZZX);
            solve_system<zz_pE>(B, U_Fq, Udenom, T_Fq, initB);
            assert_print(Udenom, ==, 1);
            U_ZZX = conv< Mat<ZZX> >(U_Fq);
        }
    }
    Unom = conv< Mat<R> >(U_ZZX);
}

template<>
inline void solve_system_local(Vec<int64_t> &B,  Mat<zz_p> &Unom, const Mat<zz_p> &T, const Vec<int64_t> &initB, const int64_t p, const ZZX f)
{
    solve_system_local_lzz_p(B, Unom, T, initB, p, f);
}
template<>
inline void solve_system_local(Vec<int64_t> &B,  Mat<ZZ_p> &Unom, const Mat<ZZ_p> &T, const Vec<int64_t> &initB, const int64_t p, const ZZX f)
{
    solve_system_local_lzz_p(B, Unom, T, initB, p, f);
}
template<>
inline void solve_system_local(Vec<int64_t> &B,  Mat<zz_pE> &Unom, const Mat<zz_pE> &T, const Vec<int64_t> &initB, const int64_t p, const ZZX f)
{
    solve_system_local_lzz_pE(B, Unom, T, initB, p, f);
}
template<>
inline void solve_system_local(Vec<int64_t> &B,  Mat<ZZ_pE> &Unom, const Mat<ZZ_pE> &T, const Vec<int64_t> &initB, const int64_t p, const ZZX f)
{
    solve_system_local_lzz_pE(B, Unom, T, initB, p, f);
}





// given two maps T : R^k -> R^n, S: R^l -> R^n, where T is injective,
// computes the columns of T not in img(T) \cap img(S)
// R = ZZ or field = zz_p or zz_pE
template<typename R>
void cokernel_intersection(Vec<int64_t> &res, const Mat<R> &T, Mat<R> &S)
{
    assert_print(T.NumRows(), ==, S.NumRows());

    Vec<int64_t> pivots, nonpivots;
    int64_t i, j;
    Mat<R> K, M;

    pivot_columns(pivots, S);
    int64_t rank_S =  pivots.length();

    //M = [ rows of T | - Base(img(S))] ^t
    M.SetDims( T.NumCols() + rank_S, T.NumRows());
    for(i = 0; i < T.NumCols(); i++)
        for(j = 0; j < T.NumRows(); j++)
            M[i][j] = T[j][i];
    for(i = 0; i < rank_S; i++)
        for(j = 0; j < T.NumRows(); j++)
            M[T.NumCols() + i][j] = -S[j][pivots[i]];

    // K = M.left_kernel()
    kernel(K, M);

    // extract the first T.NumCols() columns of K
    // this form basis for the intersection in R^n

    M.SetDims(K.NumRows(), T.NumCols());
    for(i = 0; i < K.NumRows(); i++)
        for(j = 0; j < M.NumCols(); j++)
            M[i][j] = K[i][j];

    pivot_columns(pivots, M);
    complement(res, M.NumCols(), pivots);
}

//R = zz_p or ZZ_p
template<typename R>
void cokernel_intersection_local_lzz_p(Vec<int64_t> &res, const Mat<R> &T, Mat<R> &S, const int64_t p, const ZZX f)
{
    //f = 0 or f = X
    assert_print(f == 0, or, IsX(f));
    // not an ideal implementation , but I prefer to be careful in case zz_p is being used somewhere else
    Mat<ZZ> T_ZZ = conv< Mat<ZZ> >(T);
    Mat<ZZ> S_ZZ = conv< Mat<ZZ> >(S);
    {
        zz_pPush push(p);
        Mat<zz_p> T_Fp = conv< Mat<zz_p> > (T_ZZ);
        Mat<zz_p> S_Fp = conv< Mat<zz_p> > (S_ZZ);
        cokernel_intersection<zz_p>(res, T_Fp, S_Fp);
    }
}

//R = ZZ_pE or zz_pE
template<typename R>
void cokernel_intersection_local_lzz_pE(Vec<int64_t> &res, const Mat<R> &T, Mat<R> &S, const int64_t p, const ZZX f)
{
    assert_print( f, !=, 0);
    assert_print( IsX(f), ==, false);
    Mat<ZZX> T_ZZX = conv< Mat<ZZX> >( T );
    Mat<ZZX> S_ZZX = conv< Mat<ZZX> >( S );
    {
       zz_pPush push1(p);
       {
            zz_pEPush push2(conv<zz_pX>(f));
            zz_pE Udenom;
            Mat<zz_pE> T_Fq = conv< Mat<zz_pE> >(T_ZZX);
            Mat<zz_pE> S_Fq = conv< Mat<zz_pE> >(S_ZZX);
            cokernel_intersection<zz_pE>(res, T_Fq, S_Fq);
        }
    }
}


template<typename R>
void cokernel_intersection_local(Vec<int64_t> &res, const Mat<R> &T, Mat<R> &S, const int64_t p, const ZZX f);


template<>
inline void cokernel_intersection_local(Vec<int64_t> &res, const Mat<zz_p> &T, Mat<zz_p> &S, const int64_t p, const ZZX f)
{
    cokernel_intersection_local_lzz_p(res, T, S, p, f);
}

template<>
inline void cokernel_intersection_local(Vec<int64_t> &res, const Mat<ZZ_p> &T, Mat<ZZ_p> &S, const int64_t p, const ZZX f)
{
    cokernel_intersection_local_lzz_p(res, T, S, p, f);
}

template<>
inline void cokernel_intersection_local(Vec<int64_t> &res, const Mat<zz_pE> &T, Mat<zz_pE> &S, const int64_t p, const ZZX f)
{
    cokernel_intersection_local_lzz_pE(res, T, S, p, f);
}

template<>
inline void cokernel_intersection_local(Vec<int64_t> &res, const Mat<ZZ_pE> &T, Mat<ZZ_pE> &S, const int64_t p, const ZZX f)
{
  cokernel_intersection_local_lzz_pE(res, T, S, p, f);
}



#endif

