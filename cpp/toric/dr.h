// Copyright 2017 Edgar Costa
// See LICENSE file for license details.


#ifndef DR_H
#define DR_H

#include "solve_system.h"
#include "tools.h"


#include <cstdint>
#include <assert.h>
#include <stdio.h>//needed?
#include <map>
#include <iostream>//needed?
#include <fstream>//needed?
#include <cstring>//needed?

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
template<class R>
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
        //returns true if v \in k*P  and w \notin k*P for some k, i.e., min_P(v) < min_P(w)
        bool minPcomp(const Vec<int64_t> &v, const Vec<int64_t> &w){ return min_P(v) < min_P(w); }
        
        
        // tuple_list[d] = integral points in d*P 
        // tuple_int_list[d] = integral interior points in d*P
        // we store them as (n + 1) tuples (d, w) where w \in d*P or equivalently w/d \in P
        // d \leq n + 1
        Vec< Vec<int64_t> > tuple_list;
        Vec< Vec<int64_t> > tuple_int_list;
        
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
        Vec< Vec<int64_t> > cokernels_J_basis, cokernels_I_basis, basis_dr_Y, basis_dr_X;
        // reverse maps
        Vec< map< Vec<int64_t>, int64_t, vi64less> > cokernels_J_basis_dict, cokernels_I_basis_dict, basis_dr_Y_dict, basis_dr_X_dict;


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
        Vec< Mat<int64_t> > inclusion_matrices;
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
        Mat<int64_t> proj_X, proj_notX;
        //computes  proj_X and proj_notX
        void init_proj();



        
        /*
         * Functions
         */
        
        /*
         * Constructors
         */
        dr(){};
        dr(const int64_t &p, const map< Vec<int64_t>, R, vi64less> &f, const Mat<int64_t> &AP, const Mat<int64_t> &bP,  int64_t verbose = 0, bool minimal = false){ init(p, f, AP, bP, verbose, minimal); }
        virtual ~dr(){};
        void init(const int64_t &p, const map< Vec<int64_t>, R, vi64less> &f, const Mat<int64_t> &AP, const Mat<int64_t> &bP,  int64_t verbose = 0, bool minimal = false);

        // computes matrix of the map
        // (H0, \dots Hn) ---> h0 * f + \sum_\lambda Hi * \partial_\lambda f
        // where Hi \in P_(d - 1), \lambda runs over a dual basis
        void matrix_J(Mat<R> &result, int64_t d);


        void solve_system(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<R> &T, const Vec<int64_t> &initB);


        /*
         * Test functions
         */
        void test_inclusion_matrices();
        void test_last_reduction();
    
};



template<class R>
void dr<R>::init(const int64_t &p, const map< Vec<int64_t>, R, vi64less> &f, const Mat<int64_t> &AP, const Mat<int64_t> &bP,  int64_t verbose, bool minimal)
{
    this->p = p;
    this->f = f;
    this->AP = AP;
    this->bP = bP;
    this->verbose = verbose;
    n = f.begin()->first.length();
    ring<R>(precision, fE, modulus, p);
    if( verbose > 0 )
    {
        cout<<"n = "<<n;
        cout<<" p = "<<p;
        cout<<" precision = "<<precision;
        cout<<endl;
        cout<<"AP = \n"<<AP<<endl;
        cout<<"bP = "<<bP<<endl;
        cout <<"f = \n";
        cout <<= f;
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
}

template<class R>
void dr<R>::init_tuples()
{
    Vec< Vec<int64_t>> fkeys;
    Vec<int64_t> v;
    int64_t i;
    v.SetLength(n);
    typename map< Vec<int64_t>, R, vi64less>::const_iterator fit;
    for(fit = f.begin(); fit != f.end(); fit++)
    {
        //v = fit->first[1:]
        for(i = 0; i < n; i++)
            v[i] = fit->first[i+1];
        fkeys.append(v);
    }
    
    Vec< Vec<int64_t> > local_tuple_list, local_tuple_int_list;
    integral_points(local_tuple_list, local_tuple_int_list, AP, bP, fkeys, n + 2);
    
    assert(local_tuple_list[0].length() == 1);
    assert(local_tuple_int_list[0].length() == 0 );
    //copy local data to object data
    tuple_list.SetLength(n+1);
    tuple_int_list.SetLength(n+1);
    for(i = 0; i < n + 2; i++)
    {
        tuple_list[i] = local_tuple_list[i];
        tuple_int_list[i] = local_tuple_int_list[i];
    }

    tuple_dict.SetLength(n+1);
    tuple_int_dict.SetLength(n+1);
    for(i = 0 ; i < n + 1; i++)
    {
        reverse_dict(tuple_dict[i], tuple_list[i]);
        reverse_dict(tuple_int_dict[i], tuple_int_list[i]);
    }
}


template<class R>
void dr<R>::init_f_power(int64_t N)
{
    if( f_power.length() == 0)
    {
        f_power.SetLength(1);
        Vec<int64_t> zero;
        zero.SetLength(n + 1);
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
}
template<class R>
void dr<R>::init_solve_and_cokernels()
{
    int64_t i, j;
    Vec<int64_t> ci;
    ci.SetLength(n + 1);
    for(i = 0; i < n + 1; i++)
        ci[i] = tuple_list[i].length();

    ZZX P, Q, tmp; // Q = 1 - X 
    SetCoeff(Q, 0, 1);
    SetCoeff(Q, 1, -1);
    // P = sum( ci[i] * (T ** i) * (1 - T)**(self.n - i) for i in range(self.n + 1))
    P = 0;
    for(i = 0; i < n + 1; i++)
    {
        //tmp = ci[i] * T^i
        LeftShift(tmp, conv<ZZX>(ci[i]), i);

        //tmp = ci[i] * T^i * (1-T)^(n-i)
        for(j = 0; j < n - i; j++)
            tmp *= Q;
        P+=tmp;
    }
    cokernels_J_dimensions = conv< Vec<int64_t> >(P);
    cokernels_J_dimensions.append(0);
    dim_J = sum(cokernels_J_dimensions); 

    if( verbose > 0 )
        cout << "dim J_i = "<<cokernels_J_dimensions<<endl;
    
    cokernels_J_basis.SetLength(n + 2);
    cokernels_J_basis_dict.SetLength(n + 2);
    for(i = n + 1; i >= 0; i--)
    {
        if(i == n + 1 and verbose > 0)
            cout << "Asserting that the hypersurface is non degenerate" << endl;

        Mat<R> J;
        matrix_J(J, i);
        if( verbose > 0)
            printf("Solving Jacobian relations at degree %ll (%ll x %ll)", i, J.NumRows(), J.NumCols());

        Vec<int64_t> B, initB;
        if( i <= n )
        {
            initB.SetLength(cokernels_I_dimensions[i]);
            for(j = 0; j < cokernels_I_dimensions[i]; i++)
                initB[j] = tuple_dict[i][cokernels_I_basis[j]];

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
            printf("Expected B.length() = %lld, got %lld", cokernels_J_dimensions[i], B.length());
            assert(false);
        }
        if( i > 0 )
            assert( solve_matrix[i].NumRows() == (n + 1)*tuple_list[i-1].length() + B.length() );
        
        // copy local data to object data
        cokernels_J_basis[j].SetLength(B.length());
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

    basis_dr_Y[0].SetLength[0];
    basis_dr_Y[1] = tuple_list[1];
    basis_dr_Y_dict[1] = tuple_dict[1];

    basis_dr_X[0].SetLength[0];
    basis_dr_X[1] = tuple_int_list[1];
    basis_dr_X_dict[1] = tuple_int_dict[1];

    for(i = 2; i < n + 1; i++)
    {
        basis_dr_Y[i] = cokernels_J_basis[i];
        reverse_dict(basis_dr_Y_dict[i], basis_dr_Y[i]);

        basis_dr_X[i] = cokernels_I_basis[i];
        reverse_dict(basis_dr_X_dict[i], basis_dr_X[i]);


    }

}

template<class R>
void dr<R>::init_cokernels_I_basis()
{
    cokernels_I_dimensions.SetLength(n + 1);
    cokernels_I_basis.SetLength(n + 1);
    cokernels_I_basis_dict.SetLength(n + 1);
    for(int64_t d = 1; d < n + 1; d++)
    {
        init_cokernels_I_basis(d);
        reverse_dict(cokernels_I_basis_dict[d], cokernels_I_basis[d]);
        cokernels_I_dimensions[d] = cokernels_I_basis[d].length();
    }
    dim_I = sum(cokernels_I_dimensions);
    if(verbose > 0)
        cout << "dim I_i = "<<cokernels_I_dimensions<<endl;

}

template<class R>
void dr<R>::init_cokernels_I_basis(int64_t d)
{
    assert( p != 0 );
    //I_d = (S^* + J(f) / J(f))_d
    Mat<R> JR;
    Vec<int64_t> pivots;
    int64_t i, j;
    int64_t rank_J, dim_Sintd;
    matrix_J(JR, d);
    dim_Sintd = tuple_int_list[d].length();

    //switch to a field
    if(fE == 1)
    {
        Mat<ZZ> JZZ = conv< Mat<ZZ> >(JR);
        {
            zz_pPush push(p);

            Mat<zz_p> J, T, K;
            J = conv< Mat<zz_p> >(JZZ);
            pivot_columns(pivots, J);
            rank_J = pivots.length();

            // Computing the intersection
            // (S^* \cap J(f))_d in S_d
    
            // T = [ Base(S^*_d) | - Base(img(MJ) ]^t
            T.SetDims(  dim_Sintd + rank_J, J.NumRows());
            for(i = 0; i < J.NumRows(); i++)
                for(j = 0; j < pivots.length(); j++)
                    T[rank_J + j][i] = -J[i][pivots[j]];
            for(i = 0; i < dim_Sintd; i++)
                T[i][ tuple_dict[d][tuple_int_list[i]] ] = 1;

    
            // T.left_kernel()
            kernel(K, T); 
            // extract the first dim_SintD columns of K
            // this is a basis for (S^* \cap J(f))_d in S_d
            T.SetDims(K.NumRows(), dim_Sintd);
            for(i = 0; i < K.NumRows(); i++)
                for(j = 0; j < dim_Sintd; j++)
                    T[i][j] = K[i][j];
            
            // cokernels_I_basis[d] = T.nonpivots(), the monomials not in the intersection
            pivot_columns(pivots, T);
            complement(cokernels_I_basis[d], dim_Sintd, pivots);
        }
    }
    else
    {
        //fE != 1 ==> R = zz_pE or ZZ_pE
        Mat<ZZX> JZZX = conv< Mat<ZZX> >(JR);
        {
            zz_pPush push(p);
            {
                zz_pEPush push(conv<zz_pX>(fE));

                Mat<zz_pE> J, T, K;
                J = conv< Mat<zz_pE> >(JZZX);
                pivot_columns(pivots, J);
                rank_J = pivots.length();

                // Computing the intersection
                // (S^* \cap J(f))_d in S_d
        
                // T = [ Base(S^*_d) | - Base(img(MJ) ]^t
                T.SetDims(  dim_Sintd + rank_J, J.NumRows());
                for(i = 0; i < J.NumRows(); i++)
                    for(j = 0; j < pivots.length(); j++)
                        T[rank_J + j][i] = -J[i][pivots[j]];
                for(i = 0; i < dim_Sintd; i++)
                    T[i][ tuple_dict[d][tuple_int_list[i]] ] = 1;

        
                // T.left_kernel()
                kernel(K, T); 
                // extract the first dim_SintD columns of K
                // this is a basis for (S^* \cap J(f))_d in S_d
                T.SetDims(K.NumRows(), dim_Sintd);
                for(i = 0; i < K.NumRows(); i++)
                    for(j = 0; j < dim_Sintd; j++)
                        T[i][j] = K[i][j];
                
                // cokernels_I_basis[d] = T.nonpivots(), the monomials not in the intersection
                pivot_columns(pivots, T);
                complement(cokernels_I_basis[d], dim_Sintd, pivots);

            }
        }
    }
}
//deal with the case R = ZZ
template <>
inline void dr<ZZ>::init_cokernels_I_basis(int64_t d)
{
    assert( p == 0 );
    //I_d = (S^* + J(f) / J(f))_d
    Mat<ZZ> J;
    Vec<int64_t> pivots;
    int64_t i, j;
    int64_t rank_J, dim_Sintd;
    matrix_J(J, d);
    dim_Sintd = tuple_int_list[d].length();


    Mat<ZZ> T, K;
    pivot_columns(pivots, J);
    rank_J = pivots.length();

    // Computing the intersection
    // (S^* \cap J(f))_d in S_d

    // T = [ Base(S^*_d) | - Base(img(MJ) ]^t
    T.SetDims(  dim_Sintd + rank_J, J.NumRows());
    for(i = 0; i < J.NumRows(); i++)
        for(j = 0; j < pivots.length(); j++)
            T[rank_J + j][i] = -J[i][pivots[j]];
    for(i = 0; i < dim_Sintd; i++)
        T[i][ tuple_dict[d][tuple_int_list[i]] ] = 1;


    // T.left_kernel()
    ZZ det2;
    int64_t r = image(det2, T, K);
    // m = T.NumRows()
    // the first m - r rows of K form a basis for the left kernel of T

    // extract the first  m - r rows and dim_SintD columns of K
    // this is a basis for (S^* \cap J(f))_d in S_d
    T.SetDims(T.NumRows() - r, dim_Sintd);
    for(i = 0; i < K.NumRows() - r; i++)
        for(j = 0; j < dim_Sintd; j++)
            T[i][j] = K[i][j];
    
    // cokernels_I_basis[d] = T.nonpivots(), the monomials not in the intersection
    pivot_columns(pivots, T);
    complement(cokernels_I_basis[d], dim_Sintd, pivots);

}




template<class R>
void dr<R>::matrix_J(Mat<R> &result, int64_t d)
{
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
}



template<class R>
void dr<R>::solve_system(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<R> &T, const Vec<int64_t> &initB)
{
    ::solve_system_padic<R>(B, Unom, Udenom, T, initB, p, precision, fE);
}
template <> inline void dr<ZZ>::solve_system(Vec<int64_t> &B, Mat<ZZ> &Unom, ZZ &Udenom, const Mat<ZZ> &T, const Vec<int64_t> &initB)
{
    ::solve_system(B, Unom, Udenom, T, initB);
}


/*
 * rho_i,w : P_d -> P_{d-1} 
 * pi_i : P_d --> J_d
 * m w g \omega/f^{m + 1} = w \rho_{d, w} (g) \omega/f^m + m w \pi_i(g) \omega/f^{m+1} in PH^{n-1}(Y)
 * w \in P_m and m > 0
 */

template<class R>
void dr<R>::init_rho_and_pi_matrices()
{
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
        RHO[d].SetLength(n + 2);
        // rho is represented a linear polynomial in W0, W1, ..., WN variables, that correspond to w0, ..., wn
        // RHO[0] is the constant term
        // RHO[i+1] is the Wi term, where W0 corresponds to the degree

        Mat<R> &PI = pi[d];

        R &D = solve_denom[d];
        Mat<R> Ucols;
        tranpose(Ucols, solve_matrix[d]);

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
}

template<class R>
void dr<R>::init_inclusion_matrices()
{
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
}

template<class R>
void dr<R>::test_inclusion_matrices()
{
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

    
}

template<class R>
void dr<R>::init_last_reduction()
{
    //M: P_1 + .... + P_{n + 1} --- > P_1 + J_2 + ... + J_n
    // (i - 1)! ai \in P^{i} ---> P_1 + J_2 +... + J_i
    Mat<R> &M = last_reduction;
    R &D = last_reduction_den;

    int64_t i, j, k;

    Vec<int64_t> shift_rows, shift_columns;

    shift_rows.SetLength(n + 1);
    shift_columns.SetLength(n + 1);
    // shift_rows[k] = dim(P_1 + J_2 + ... + J_{k - 1})
    // shift_columns[k] = dim(P_1 + P_2 + .... + P_{k - 1}
    int64_t total_columns, total_rows;
    total_columns = 0;
    total_rows = 0;
    for(k = 1; k < n + 1; k++)
    {
        shift_columns[k] = total_columns;
        shift_rows[k] = total_rows;
        total_columns += tuple_list[k].length();
        total_rows += basis_dr_Y[k].length();
    }
    M.SetDims(total_rows, total_columns);

    //P_1 --> P_1 is the identity map
    for(k = 1; k < tuple_list[1].length(); k++)
        M[k][k] = 1;
    
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
            G.SetLength(total_rows);
            for(i = 0; i < rho[k][0].NumRows(); i++)
                G[i + shift_columns[k - 1]] = rho[k][0][i][j];
            
            H = (M * G)/ rho_den[k];

            for(i = 0; i < H.length(); i++)
                M[i][ shift_columns[k] + j] += H[i];
        }
        factorial *= k;
    }
}

template<class R>
void dr<R>::test_last_reduction()
{
    if( verbose > 0 )
        cout << "dr.test_last_reduction()" << endl;

    Vec< Vec<int64_t> > B;
    int64_t i, j, k;
    
    for(i = 0; i < basis_dr_Y.length(); i++)
        B.append( basis_dr_Y[i] );

    for(i = 0; i < B.length(); i++)
    {
        Vec<int64_t> &v = B[i];
        Vec<R> expected;
        expected.SetLength(B.length());
        expected[i] = last_reduction_den;
        for(k = 0; k < n + 2 - v[0]; k++)
        {
            Vec<R> G;
            G.SetLength(last_reduction.NumCols());
                //FIXME
        }


    }

    /*

        for l, v in enumerate(B):
            expected = vector([0] * len(B));
            expected[l] = self.last_reduction[1] ;
            for k in range(self.n + 2 - v[0]):
                G = vector([0] * sum( map(len, self.tuple_list[1:self.n + 2])))
                shift = sum( map(len, self.tuple_list[1:k + v[0]]));
                for mon, cf in self.get_f_power(k).iteritems():
                    G[shift + self.tuple_dict[k + v[0]][immutable(mon + v)]] += cf;

                assert self.R(factorial(k + v[0] - 1)) * expected ==  (self.last_reduction[0] * G)
    */
}


template<class R>
void dr<R>::init_proj()
{
    //FIXME
    /*
     * ef project_dR(self):
        if self.verbose > 2:
            print "dr.project_dR()"

        FdR = Matrix(ZZ, self.dim_dR, self.dim_dR_T)
        FnotdR = Matrix(ZZ, self.dim_dR_T, self.dim_dR_T)
        B = [];
        for TL in self.basis_dR_T:
            B += TL
        
        for k, ck in enumerate(B):
            if ck in self.basis_dR[ck[0]]:
                FdR[sum(map(len,  self.basis_dR[:ck[0]])) + self.basis_dR_dict[ck[0]][ck], k] = 1
            else:
                FnotdR[k, k] = 1
        return FdR, FnotdR

        */
}
#endif
