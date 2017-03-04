// Copyright 2017 Edgar Costa
// See LICENSE file for license details.


#ifndef DR_H
#define DR_H

#include "linear_algebra.h"
#include "tools.h"


#include <cstdint>
#include <assert.h>
#include <stdio.h>//needed?
#include <algorithm> //sort
#include <map>
#include <iostream>//needed?
#include <fstream>//needed?
#include <cstring>//needed?
#include <sstream>


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

//R = ZZ, ZZ_p, zz_p. ZZ_pE, or zz_pe
template<typename R>
class dr{
    public:
        /*
         * Notation:
         * P a polyhedron in ZZ^n
         * f a polynomial x1, ..., xn
         * and convexhull( supp(f) U 0) =  P
         * P_k := k * P
         *
         * P^* = int(P)
         *
         * S (resp S^*), the graded ring associated to P (resp P^*)
         * S := \sum P_k 
         * S^* := \sum P^* _k
         *
         * P the toric variety associated to P
         * T (= C^* ^n) the torus in P 
         * X < P, the closure of V(f) in P
         * Y = X \cap T, hypersurface defined by f in the Torus, i.e, V(f)
         * 
         * J(f) = < f, \partial_lambda f : \lambda \in (ZZ^n)^>
         *   = < f, x_i \partial_i f : i = 1, ... n >
         * For simplicity, write F0 = f and , Fi = \partial_i f
         *
         * write J_d = (S / J(f))_d
         * write I_d = (S^* + J(f) / J(f))_d
         * 
         * then:
         * PH^{n-1} (Y) = P_1 + J_2 + ... + J_n
         * PH^{n-1} (X) = P^*_1 + I_2 + ... + I_n 
         */

        /* 
         * Attributes
         */
        map< Vec<int64_t>, R, vi64less> f; 
        int64_t n;

        /*
         * verbose = 0 --> silent
         * verbose = 1 --> minimal
         * verbose = 2 --> human readable
         * verbose > 2 --> everything
         */
        int64_t verbose;
        bool minimal;
        // Specifications that must match with the modulus that defines R
        int64_t p; // the characteristic where f comes from, = 0 if R is ZZ, else the characteristic of FFq
        int64_t precision;
        ZZX fE; // 1 if p == q else FFq.defining_polynomial()
        ZZ modulus; // p^precision
        

        // Half space representation of P
        // w \in k*P <=>  AP * v + k *bP >= 0
        Mat<int64_t> AP;
        Vec<int64_t> bP;
        // returns the minimal k such that v \in k*P
        int64_t min_P(const Vec<int64_t> &v){ return ::min_P(AP, bP, v); }
        //returns the minimal k such that v \in k*int(P)
        int64_t min_intP(const Vec<int64_t> &v){ return ::min_intP(AP, bP, v); }
                


        
        
        // tuple_list[d] = integral points in d*P 
        // tuple_int_list[d] = integral interior points in d*P
        // we store them as (n + 1) tuples (d, w) where w \in d*P or equivalently w/d \in P
        // d \leq n + 1
        Vec< Vec< Vec<int64_t> > > tuple_list;
        Vec< Vec< Vec<int64_t> > > tuple_int_list;
        
        // reverse maps
        Vec< map< Vec<int64_t>, int64_t, vi64less> > tuple_dict;
        Vec< map< Vec<int64_t>, int64_t, vi64less> > tuple_int_dict;
        
        // computes all the above
        void init_tuples();
        
        // f_power[i] = f^i
        Vec<  map< Vec<int64_t>, R, vi64less> > f_power;
        // if needed computes f^N and adds it to the f_power list
        void init_f_power(int64_t N);
    
        // solve_matrix[d] represents the map
        // P_d ---> P_{d-1}^n + J_d 
        //   g ---> (gi, ci)
        // such that
        // g = \sum gi * Fi + cokernels 
        Vec< Mat<R> > solve_matrix;
        Vec< R > solve_denom;
        
        // coKernels_J_dimensions[i] = \dim J_i
        Vec<int64_t> cokernels_J_dimensions, cokernels_I_dimensions;
        int64_t dim_J, dim_I; //sum(coKernels_*_dimensions)
        
        // cokernels_*_basis[i] = basis for *_i
        // constructed such that I_i < J_i
        // basis_dr_Y[i] = J_i or P_1 if i == 1
        // basis_dr_X[i] = I_i or P^*_1 if i == 1 FIXME: P^*_1 = I_1
         Vec< Vec< Vec<int64_t> > >cokernels_J_basis, cokernels_I_basis, basis_dr_Y, basis_dr_X;
        // reverse maps
        Vec< map< Vec<int64_t>, int64_t, vi64less> > cokernels_J_basis_dict, cokernels_I_basis_dict, basis_dr_Y_dict, basis_dr_X_dict;
        int64_t dim_dr_Y, dim_dr_X;


        //  computes cokernels_I_basis*
        void init_cokernels_I_basis();
        void init_cokernels_I_basis(int64_t d);

        // computes solve_*, cokernels_J_basis*, basis_dr_*
        // must be called after  init_cokernels_I_basis()
        void init_solve_and_cokernels();
        
        



        // rho_i,w : P_i --> P_{i-1}
        // m w g \omega/f^{m +1} = w \rho_{i, w} (g) \omega/f^m + m w \pi_i(g) \omega/f^{m+1}
        // in PH^{n-1}(Y)
        //
        // rho[i] stores rho_i as one degree polynomial in (W_1, ..., W_n) with matrix coefficients
        // rho[i][j] represents the coefficient of rho_i associated to W_{j+1}
        // e.g. rho[i][0] is the constant coefficient
        Vec< Vec< Mat<R> > > rho;
        Vec< R > rho_den;

        // pi[i] represents the matrix
        // pi_i: P_i --> J(f)_{i-1}
        Vec< Mat<R> > pi;
        Vec< R > pi_den;

        // computes rho[d] and pi[d]
        void init_rho_and_pi_matrices();

        // inclusion_matrices[u] is the matrix that represents:
        // I_u : J_0 + ... + J_n --> P_1 + P_2 + ... + P_{n+1}
        // where ai -- >  u * a_i \in P_{i + 1} and u \in P1
        Vec< Mat<R> > inclusion_matrices;
        //compute inclusion_matrices
        void init_inclusion_matrices();

        // last_reduction represents
        // M: P_1 + .... + P_{n + 1} --- > P_1 + J_2 + ... + J_n
        // (i - 1)! ai \in P^{i} ---> P_1 + J_2 +... + J_i
        Mat<R> last_reduction;
        R last_reduction_den;
        //computes the last_reduction
        void init_last_reduction();


        // projection of matrices 
        // proj_X :  PH^{n-1} (Y) ---> PH^{n-1} (X)
        // proj_notX :  PH^{n-1} (Y) ---> PH^{n-1} (X)^{perp}
        Mat<R> proj_X, proj_notX;
        //computes  proj_X and proj_notX
        void init_proj();



        
        /*
         * Functions
         */
        
        /*
         * Constructors
         */
        dr(){};
        dr(const char* input, const int64_t &verbose = 0, const bool &minimal = false);
        dr(const int64_t &p, const map< Vec<int64_t>, R, vi64less> &f, const Mat<int64_t> &AP, const Vec<int64_t> &bP,  const int64_t &verbose = 0, const bool &minimal = false){ init(p, f, AP, bP, verbose, minimal); }
        virtual ~dr(){};
        void init(const int64_t &p, const map< Vec<int64_t>, R, vi64less> &f, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const  int64_t &verbose = 0, const bool &minimal = false);

        // computes matrix of the map
        // (H0, \dots Hn) ---> h0 * f + \sum_\lambda Hi * \partial_\lambda f
        // where Hi \in P_(d - 1), \lambda runs over a dual basis
        void matrix_J(Mat<R> &result, int64_t d);


        void solve_system(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<R> &T, const Vec<int64_t> &initB);
        void cokernel_intersection(Vec<int64_t>  &res, Mat<R> &T, Mat<R> &S);


        /*
         * Test functions
         */
        void test_inclusion_matrices();
        void test_last_reduction();
        void test_all();
    
};

template<typename R>
dr<R>::dr(const char* input, const int64_t &verbose, const bool &minimal)
{
    /*
     input: 
     p
     f.keys()
     f.values()
     AP
     bP

     for example:
     '17 \n[[0 2][0 0][3 0][2 0][0 1]]\n[1 12 16 1 1] \n[[ 0  1][ 1  0][-2 -3]] \n[0 0 6] \n'
     */

    stringstream buffer;
    buffer << input;
    int64_t local_p;
    map< Vec<int64_t>, R, vi64less> local_f;
    Mat<int64_t> local_AP;
    Vec<int64_t> local_bP;
    buffer >> local_p;
    buffer >> local_f;
    buffer >> local_AP;
    buffer >> local_bP;
    init(local_p, local_f, local_AP, local_bP, verbose, minimal);
}


template<typename R>
void dr<R>::init(const int64_t &p, const map< Vec<int64_t>, R, vi64less> &f, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const int64_t &verbose, const bool &minimal)
{
    if(verbose > 2)
        cout<<"dr::init()"<<endl;
    this->p = p;
    n = f.begin()->first.length();
    typename map< Vec<int64_t>, R, vi64less>::const_iterator fit;
    Vec<int64_t> v;
    v.SetLength(n + 1);
    v[0] = 1;
    for(fit = f.begin(); fit != f.end(); fit++)
    {
        for(int64_t i = 0; i < n; i++)
            v[i + 1] = fit->first[i];
        this->f[v] = fit->second;
    }

    this->AP = AP;
    this->bP = bP;
    this->verbose = verbose;
    ring<R>(precision, fE, modulus, p);
    if( verbose > 0 )
    {
        cout<<"n = " << this->n;
        cout<<" p = " << this->p;
        cout<<" precision = " << this->precision;
        cout<<endl;
        cout<<"AP = \n" << this->AP<<endl;
        cout<<"bP = " << this->bP<<endl;
        cout << "f = \n";
        cout <<= this->f;
        cout <<endl;
    }

    init_tuples();
    init_cokernels_I_basis();
    if(not minimal)
    {
        init_solve_and_cokernels();
        init_rho_and_pi_matrices();
        init_last_reduction();
        init_proj();
    }
    if(verbose > 2)
        cout<<"dr::init() done"<<endl;
}

template<typename R>
void dr<R>::init_tuples()
{
    if(verbose > 2)
        cout<<"dr::init_tuples() begin"<<endl;

    Vec< Vec<int64_t>> fkeys;
    fkeys.SetLength(0);
    Vec<int64_t> v;
    int64_t i;
    v.SetLength(n);
    v[0] = 1;
    typename map< Vec<int64_t>, R, vi64less>::const_iterator fit;
    for(fit = f.begin(); fit != f.end(); fit++)
    {
        //v = fit->first[1:]
        for(i = 0; i < n; i++)
            v[i] = fit->first[i + 1];
        fkeys.append(v);
    }
    integral_points(tuple_list, tuple_int_list, AP, bP, fkeys, n + 2);
    
    assert_print(tuple_list[0].length(), ==, 1);
    assert_print(tuple_int_list[0].length(), ==, 0 );
    
    tuple_dict.SetLength(n + 2);
    tuple_int_dict.SetLength(n + 2);
    
    //returns true if v \in k*P  and w \notin k*P for some k, i.e., min_P(v) < min_P(w)
    auto minPless = [=](const Vec<int64_t> &v, const Vec<int64_t> &w){return ::min_intP(AP, bP, v) < ::min_intP(AP, bP, w);};
 


    for(i = 0 ; i < n + 2; i++)
    {
        sort(tuple_int_list[i].begin(), tuple_int_list[i].end(), vi64less() );
        stable_sort(tuple_int_list[i].begin(), tuple_int_list[i].end(), minPless);
        stable_sort(tuple_list[i].begin(), tuple_list[i].end(), minPless);
        reverse_dict(tuple_dict[i], tuple_list[i]);
        reverse_dict(tuple_int_dict[i], tuple_int_list[i]);
    }
    if(verbose > 2)
        cout<<"dr::init_tuples() end"<<endl;
}


template<typename R>
void dr<R>::init_f_power(int64_t N)
{
    if(verbose > 2)
        cout<<"dr::init_f_power("<<N<<")"<<endl;
    if( f_power.length() == 0)
    {
        f_power.SetLength(1);
        Vec<int64_t> zero;
        zero.SetLength(n + 1, 0);
        f_power[0][zero] = R(1);
    }
    
    if( f_power.length() < N + 1)
    {
        int64_t k;
        k = f_power.length();
        f_power.SetLength(N + 1);
        typename map< Vec<int64_t>, R, vi64less>::const_iterator fit, git;
        Vec<int64_t> u;

        for(; k < N + 1; k++)
        {
             for(git = f_power[k - 1].begin(); git != f_power[k - 1].end(); git++)
             {
                for(fit = f.begin(); fit != f.end(); fit++)
                {
                    u = fit->first + git->first;
                    f_power[k][u] += fit->second * git->second;
                }
             }

        }
    }
    if(verbose > 2)
        cout<<"dr::init_f_power("<<N<<") done"<<endl;
}
template<typename R>
void dr<R>::init_solve_and_cokernels()
{
    if(verbose > 2)
        cout<<"dr::init_solve_and_cokernels()"<<endl;
    int64_t i, j;
    Vec<ZZ> ci;
    // ci = rank(P_i)
    ci.SetLength(n + 2);
    for(i = 0; i < n + 2; i++)
        ci[i] = tuple_list[i].length();
    
    // \sum rank(J_i) T^i = (1- T)^(n+1) \sum rank(P_i) T^i
    ZZX P, Q, tmp; // Q = 1 - T
    SetCoeff(Q, 0, 1);
    SetCoeff(Q, 1, -1);
    // P = (1- T)^(n+1) \sum rank(P_i) T^i
    P = conv<ZZX>(ci);
    for(i = 0; i < n + 1; i++)
        P *= Q;
    cokernels_J_dimensions.SetLength(n + 2, 0);
    assert_print(coeff(P, n + 1), ==, 0);
    for(i = 0; i < n + 2; i++)
        cokernels_J_dimensions[i] = conv< long >( coeff(P, i) );
    dim_J = sum(cokernels_J_dimensions); 

    if( verbose > 1 )
        cout << "dim J_i = "<<cokernels_J_dimensions<<endl;
    
    cokernels_J_basis.SetLength(n + 2);
    cokernels_J_basis_dict.SetLength(n + 2);
    solve_matrix.SetLength(n + 2);
    solve_denom.SetLength(n + 2);
    for(i = n + 1; i >= 0; i--)
    {
        if(i == n + 1 and verbose > 1)
            cout << "Asserting that the hypersurface is non degenerate" << endl;

        Mat<R> J;
        matrix_J(J, i);
        if( verbose > 1)
            printf("Solving Jacobian relations at degree %ld (%ld x %ld)\n", (long)i, (long)J.NumRows(), (long)J.NumCols());

        Vec<int64_t> B, initB;
        if( i <= n )
        {
            initB.SetLength(cokernels_I_dimensions[i]);
            for(j = 0; j < cokernels_I_dimensions[i]; j++)
                initB[j] = tuple_dict[i][ cokernels_I_basis[i][j] ];

        }
        else
        {
            initB.SetLength(0);
        }
        solve_system(B, solve_matrix[i], solve_denom[i], J, initB);
        if(B.length() != cokernels_J_dimensions[i])
        {
            if(i == n + 1)
                cout<<"The hypersurface is not nondegenerate"<<endl;
            printf("Expected B.length() = %ld, got %ld", (long) cokernels_J_dimensions[i], (long) B.length());
            assert(false);
        }
        if( i > 0 )
            assert( solve_matrix[i].NumRows() == (n + 1)*tuple_list[i-1].length() + B.length() );
        
        
        // copy local data to object data
        cokernels_J_basis[i].SetLength(B.length());
        for(j = 0; j < B.length(); j++)
            cokernels_J_basis[i][j] = tuple_list[i][B[j]];

        reverse_dict(cokernels_J_basis_dict[i], cokernels_J_basis[i]);
    }

    // PH^{n-1} (Y) = P_1 + J_2 + ... + J_n
    // PH^{n-1} (X) = P^*_1 + I_2 + ... + I_n 
    basis_dr_Y.SetLength(n+1);
    basis_dr_Y_dict.SetLength(n+1);

    basis_dr_X.SetLength(n+1);
    basis_dr_X_dict.SetLength(n+1);

    basis_dr_Y[0].SetLength(0);
    basis_dr_Y[1] = tuple_list[1];
    basis_dr_Y_dict[1] = tuple_dict[1];
    dim_dr_Y = tuple_list[1].length();

    basis_dr_X[0].SetLength(0);
    basis_dr_X[1] = tuple_int_list[1];
    basis_dr_X_dict[1] = tuple_int_dict[1];
    dim_dr_X =  tuple_int_list[1].length();


    for(i = 2; i < n + 1; i++)
    {
        basis_dr_Y[i] = cokernels_J_basis[i];
        reverse_dict(basis_dr_Y_dict[i], basis_dr_Y[i]);
        dim_dr_Y += cokernels_J_basis[i].length();


        basis_dr_X[i] = cokernels_I_basis[i];
        reverse_dict(basis_dr_X_dict[i], basis_dr_X[i]);
        dim_dr_X += cokernels_I_basis[i].length();
    }
    if(verbose > 2)
        cout<<"dr::init_solve_and_cokernels() done"<<endl;
}

template<typename R>
void dr<R>::init_cokernels_I_basis()
{
    if(verbose > 1)
        cout<<"dr::init_cokernels_I_basis()"<<endl;
    cokernels_I_dimensions.SetLength(n + 1, 0);
    cokernels_I_basis.SetLength(n + 1);
    cokernels_I_basis_dict.SetLength(n + 1);
    for(int64_t d = 1; d < n + 1; d++)
    {
        init_cokernels_I_basis(d);
        reverse_dict(cokernels_I_basis_dict[d], cokernels_I_basis[d]);
        cokernels_I_dimensions[d] = cokernels_I_basis[d].length();
    }
    dim_I = sum(cokernels_I_dimensions);
    if(verbose > 1)
        cout << "dim I_i = "<<cokernels_I_dimensions<<endl;
    if(verbose > 2)
        cout<<"dr::init_cokernels_I_basis() done"<<endl;


}

template<typename R>
void dr<R>::init_cokernels_I_basis(int64_t d)
{
    if(verbose > 2)
        cout<<"dr::init_cokernels_I_basis("<<d<<")"<<endl;

    //I_d = (S^* + J(f) / J(f))_d
    Mat<R> J;
    Vec<int64_t> nonpivots;
    int64_t i, dim_Sintd;
    matrix_J(J, d);
    dim_Sintd = tuple_int_list[d].length();

    //Matrinx representing S^*_d in S_d
    Mat<R> T;
    T.SetDims(J.NumRows(), dim_Sintd);
    for(i = 0; i < dim_Sintd; i++)
        T[tuple_dict[d][tuple_int_list[d][i]]][i] = 1;

    //computes the cokernel of the map (img(J) \cap S^* )_d -> S^*_d
    cokernel_intersection(nonpivots, T, J);
    cokernels_I_basis[d].SetLength(nonpivots.length());
    for(i = 0; i < nonpivots.length(); i++)
        cokernels_I_basis[d][i] = tuple_int_list[d][nonpivots[i]];

    if(verbose > 2)
        cout<<"dr::init_cokernels_I_basis("<<d<<") done"<<endl;
}




template<typename R>
void dr<R>::matrix_J(Mat<R> &result, int64_t d)
{
    if(verbose > 2)
        cout<<"dr:::matrix_J(-, "<<d<<")"<<endl;


    //free storage and make 0 x 0 
    result.kill();
    assert( d >= 0 );
    if(d == 0)
    {
        result.SetDims( tuple_list[d].length(),0);
        return;
    }
    result.SetDims( tuple_list[d].length(), (n + 1) * tuple_list[d - 1].length());
    
    typename map< Vec<int64_t>, R, vi64less>::const_iterator fit;
    for(int64_t j = 0; j <  tuple_list[d - 1].length(); j++)
    {
        Vec<int64_t> v = tuple_list[d - 1][j];
        for(fit = f.begin(); fit != f.end(); fit++)
            for( int64_t i = 0; i < n + 1; i++)
                result[tuple_dict[d][fit->first + v]][ i * tuple_list[d - 1].length() + j] = fit->first[i] * fit->second;
    }
    if(verbose > 2)
        cout<<"dr:::matrix_J(-, "<<d<<") done"<<endl;
}



template<typename R>
void dr<R>::solve_system(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<R> &T, const Vec<int64_t> &initB)
{
    ::solve_system_padic<R>(B, Unom, Udenom, T, initB, p, precision, fE);
}
template <> inline void dr<ZZ>::solve_system(Vec<int64_t> &B, Mat<ZZ> &Unom, ZZ &Udenom, const Mat<ZZ> &T, const Vec<int64_t> &initB)
{
    ::solve_system(B, Unom, Udenom, T, initB);
}

template<typename R>
void dr<R>::cokernel_intersection(Vec<int64_t>  &res, Mat<R> &T, Mat<R> &S)
{
    ::cokernel_intersection_local(res, T, S, p, fE);
}
template<>
inline void dr<ZZ>::cokernel_intersection(Vec<int64_t>  &res, Mat<ZZ> &T, Mat<ZZ> &S)
{
    ::cokernel_intersection(res, T, S);
}


/*
 * rho_i,w : P_d -> P_{d-1} 
 * pi_i : P_d --> J_d
 * m w g \omega/f^{m + 1} = w \rho_{d, w} (g) \omega/f^m + m w \pi_i(g) \omega/f^{m+1} in PH^{n-1}(Y)
 * w \in P_m and m > 0
 */

template<typename R>
void dr<R>::init_rho_and_pi_matrices()
{
    if(verbose > 2)
        cout<<"dr:init_rho_and_pi_matrices()"<<endl;


    int64_t i, j, k, d;
    rho.SetLength(n + 2);
    rho_den.SetLength(n + 2);
    pi.SetLength(n + 2);
    pi_den.SetLength(n + 2);

    //deal with d = 0
    // rho[0] = 0 map, the codomain doesn't exist...
    // pi[0] is the identity on {0}
    rho.SetLength(n+2);
    rho_den[0] = 1;
    pi[0].SetDims(1, 1);
    pi[0][0][0] = 1;
    pi_den[0] = 1;
    

    for(d = 1; d < n + 2; d++)
    {
        Vec<Mat<R>> &RHO = rho[d];
        RHO.SetLength(n + 2);
        // rho is represented a linear polynomial in W0, W1, ..., WN variables, that correspond to w0, ..., wn
        // RHO[0] is the constant term
        // RHO[i+1] is the Wi term, where W0 corresponds to the degree

        Mat<R> &PI = pi[d];

        Mat<R> Ucols;
        transpose(Ucols, solve_matrix[d]);

        pi_den[d] = solve_denom[d];
        rho_den[d] = solve_denom[d];

        Vec< Vec<int64_t> > &Gtuples = tuple_list[d];
        Vec< Vec<int64_t> > &Htuples = tuple_list[d - 1];
        int64_t len_Htuples = Htuples.length();

        for(i = 0; i < n + 2; i++)
            RHO[i].SetDims( Htuples.length(), Gtuples.length() );

        PI.SetDims( cokernels_J_basis[d].length(), Gtuples.length() );

        for(j = 0; j < Gtuples.length(); j++)
        {
            /*
             * Let v = Gtuples[j]
             * We have
             * v =  \sum_i Hi \partial_i f + cokernels
             * where \partial_0 f = f
             * and \partial_i f = df / dxi
             */

            Vec<R> &H = Ucols[j];

            // deal with cokernels
            for(i = 0; i < cokernels_J_basis[d].length(); i++)
                PI[i][j] = Ucols[j][ (n + 1) * len_Htuples + i ];
        
            /*
             * no need to do anything special while dealing with the term:        
             * w  H0  f
             *
             * one one hand: 
             *      m w H0 * f \omega/ f^{m + 1} = m w H0 \omega/f^{m}
             * On the other hand 
             *      w H0 f \omega/ f^{m + 1} = w H0 (\partial_0 f) \omega/ f^{m + 1}
             * where \partial_0 =  \deg 
             * and
             *      \partial_0 (w H0) = (w0 + d - 1) w H0, where w0 + d = m + 1
             */

            for(i = 0; i < n + 1; i++)
            {
                /*
                 * dealing with the term
                 * w Hi * \partial_i f 
                 *  --> \partial_i (w Hi)
                 * = w ( wi Hi + \partial_i Hi)
                 */
                for(k = 0; k < len_Htuples; k++)
                {
                    // the constant term
                    RHO[0][k][j] += Htuples[k][i] * H[i*len_Htuples + k];
                    //the Wi term
                    RHO[i + 1][k][j] +=  H[i*len_Htuples + k];
                }
            }
        }
    }
    if(verbose > 2)
        cout<<"dr:init_rho_and_pi_matrices() done"<<endl;

}

template<typename R>
void dr<R>::init_inclusion_matrices()
{
    if(verbose > 2)
        cout<<"dr:init_inclusion_matrices()"<<endl;

    // for u in P1 the matrix representing
    // I_u : J_0 + ... + J_n --> P_1 + P_2 + ... + P_{n+1}
    // where J_i --> P_{i+1} corresponds to 
    //       ai -- >  u * a_i \in P_{i + 1}

    inclusion_matrices.SetLength(tuple_list[1].length());

    Vec<int64_t> shift_rows, shift_columns;
    shift_rows.SetLength(n + 1);
    shift_columns.SetLength(n + 1);
    int64_t total_columns, total_rows;
    total_columns = 0;
    total_rows = 0;
    for(int64_t k = 0; k < n + 1; k++)
    {
        shift_columns[k] = total_columns;
        shift_rows[k] = total_rows;
        total_columns += cokernels_J_dimensions[k];
        total_rows += tuple_list[k + 1].length();
    } 

    for(int64_t i = 0; i < tuple_list[1].length(); i++)
    {
        Mat<R> &M = inclusion_matrices[i];
        M.SetDims(total_rows, total_columns);

        Vec<int64_t> &u = tuple_list[1][i];
        for(int64_t k = 0; k < n + 1; k++)
        {
            
            for(int64_t j=0; j < cokernels_J_dimensions[k]; j++)
                M[shift_rows[k] + tuple_dict[k + 1][ u + cokernels_J_basis[k][j] ]][shift_columns[k] + j] = 1;
        }

    }
    if(verbose > 2)
        cout<<"dr:init_inclusion_matrices() done"<<endl;
}

template<typename R>
void dr<R>::test_inclusion_matrices()
{
    if(verbose > 2)
        cout<<"dr:test_inclusion_matrices()"<<endl;
    int64_t i, j, k;
    Vec< Vec<int64_t> > B, BJ;
    for(i = 1; i < n + 2; i++)
        B.append(tuple_list[i]);
    for(i = 0; i < n + 1; i++)
        BJ.append(cokernels_J_basis[i]);


    for(k = 0; k < tuple_list[1].length(); k++)
    {
        Vec<int64_t> &u = tuple_list[1][k];
        Mat<R> &M = inclusion_matrices[k];
        assert(M.NumCols() == BJ.length() and M.NumRows() == B.length());
        for(j = 0; j < BJ.length(); j++)
        {
            Vec<int64_t> &wj = BJ[j];
            for(i = 0; i < B.length(); i++)
            {
                if( M[i][j] != 0 )
                {
                    assert(u + wj == B[i]);
                    break;
                }

            }
        }

    }
    if(verbose > 2)
        cout<<"dr:test_inclusion_matrices() done"<<endl;
}

template<typename R>
void dr<R>::init_last_reduction()
{
    if(verbose > 2)
        cout<<"dr:init_last_reduction()"<<endl;

    //M: P_1 + .... + P_{n + 1} --- > P_1 + J_2 + ... + J_n
    // (i - 1)! ai \in P^{i} ---> P_1 + J_2 +... + J_i
    Mat<R> &M = last_reduction;
    R &D = last_reduction_den;

    int64_t i, j, k;

    Vec<int64_t> shift_rows, shift_columns;

    shift_rows.SetLength(n + 2, 0);
    shift_columns.SetLength(n + 2, 0);
    // shift_rows[k] = dim(P_1 + J_2 + ... + J_{k - 1})
    // shift_columns[k] = dim(P_1 + P_2 + .... + P_{k - 1}
    int64_t total_columns, total_rows;
    total_columns = 0;
    total_rows = 0;
    for(k = 1; k < n + 2; k++)
    {
        shift_columns[k] = total_columns;
        shift_rows[k] = total_rows;
        total_columns += tuple_list[k].length();
        if(k < n + 1)
            total_rows += basis_dr_Y[k].length();
    }
    M.SetDims(total_rows, total_columns);

    //P_1 --> P_1 is the identity map
    for(k = 0; k < tuple_list[1].length(); k++)
        M[k][k] = 1;
    D = 1;
    for(k = 2; k < n + 2; k++)
        D *= rho_den[k];

    int64_t factorial = 1;
    for(k = 2; k < n + 2; k++)
    {

        // P_k ---> P_1 + J_2 + ... + J_k
        // we compute P_k --> J_k, and P_k --> P_{k-1} with pi and rho respectively
        // then reuse M to compute P_{k-1} --> P_1 + J_2 + ... + J_{k-1}
        for(j = 0; j < tuple_list[k].length(); j++)
        {
            R c = D / pi_den[k];
            c *= factorial; // factorial = (k-1)!
            for(i = 0; i < pi[k].NumRows(); i++)
                M[shift_rows[k] + i][shift_columns[k] + j] += c * pi[k][i][j];

            // using rho_k,0 ( ej ) \in P_{k - 1}
            // rho[k][0].column(j).list() \in P_{k-1}
            Vec<R> G, H;
            // G = [0]*dim(P_1 + ... + P_{k-2}) +  rho[k][0].column(j) + [0]*dim(P_k + ... P_{n+1} 
            G.SetLength(total_columns, R(0));
            for(i = 0; i < rho[k][0].NumRows(); i++)
                G[i + shift_columns[k - 1]] = rho[k][0][i][j];
            H = (M * G);
            if( rho_den[k] != 1)
                for(i = 0; i < H.length(); i++)
                    H[i] /= rho_den[k];
            for(i = 0; i < H.length(); i++)
                M[i][ shift_columns[k] + j] += H[i];
        }
        factorial *= k;
    }
    if(verbose > 2)
        cout<<"dr:init_last_reduction() done"<<endl;
}

template<typename R>
void dr<R>::test_last_reduction()
{
    if( verbose > 2 )
        cout << "dr::test_last_reduction()" << endl;

    Vec< Vec<int64_t> > B;
    int64_t i, j, k, shift;
    
    for(i = 0; i < basis_dr_Y.length(); i++)
        B.append( basis_dr_Y[i] );

    init_f_power(n + 1);
    for(i = 0; i < B.length(); i++)
    {
        Vec<int64_t> &v = B[i];
        Vec<R> expected;
        expected.SetLength(B.length(), R(0));
        expected[i] = last_reduction_den;
        for(k = 0; k < n + 2 - v[0]; k++)
        {
            Vec<R> G;
            G.SetLength(last_reduction.NumCols(), R(0));
            shift = 0;
            for(j = 1; j < k + v[0]; j++)
                shift += tuple_list[j].length();

            typename map< Vec<int64_t>, R, vi64less>::const_iterator fkit;
            for(fkit = f_power[k].cbegin(); fkit != f_power[k].cend(); fkit++)
                G[shift + tuple_dict[k + v[0]][fkit->first + v]] += fkit->second;
            print(last_reduction); 
            assert_print(last_reduction_den * factorial<R>(k + v[0] - 1) * expected, ==, last_reduction * G);
        }
    }
    if( verbose > 2 )
        cout << "dr::test_last_reduction() done" << endl;
}


template<typename R>
void dr<R>::init_proj()
{
    if( verbose > 2)
        cout << "dr::project_dR()" << endl;

    proj_X.SetDims( dim_dr_X, dim_dr_Y);
    proj_notX.SetDims( dim_dr_Y, dim_dr_Y);

    int64_t i, j, shiftX, shiftY;
    shiftX = 0;
    shiftY = 0;
    for(i = 1; i < n + 1; i++)
    {
        for(j = 0; j < basis_dr_Y[i].length(); j++)
        {
            Vec<int64_t> &v = basis_dr_Y[i][j];
            if( basis_dr_X_dict[i].count(v) == 1)
                proj_X[shiftX + basis_dr_X_dict[i][v]][ shiftY + j] = 1;
            else
                proj_notX[shiftY + j][shiftY + j] = 1;
        }
        shiftX += basis_dr_X[i].length();
        shiftY += basis_dr_Y[i].length();

    }
    if( verbose > 2)
        cout << "dr::project_dR() done" << endl;


}
template<typename R>
void dr<R>::test_all()
{
    if( verbose > 0)
        cout << "dr::test_all()" << endl;
    test_last_reduction();
    test_inclusion_matrices();
    //FIXME
    //test_monomial_to_basis();

    if( verbose > 0)
        cout << "dr::test_all() done" << endl;
}

#endif
