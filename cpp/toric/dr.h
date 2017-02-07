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
         * f \in P
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
        int64_t characteristic;
        ZZ modulus;
        

        // Half space representation of P
        // w \in k*P <=>  AP * v + k *bP >= 0
        Vec<int64_t> AP;
        Mat<int64_t> bP;
        
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
        
        // represents the basis of the dual lattice associated to the polytope (of the polyhedron)
        // dual_basis[i] = e[i]  i < n + 1
        Vec< Vec<int64_t> > dual_basis;

        // f_power[i] = f^i
        Vec<  map< Vec<int64_t>, R, vi64less> > f_power;
        // if needed computes f^N and adds it to the f_power list
        void compute_f_power(int64_t N);
    
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
        // and basis_dr_Y[i] = J_i or P_1 if i == 1
        Vec< Vec<int64_t> > cokernels_J_basis, cokernels_I_basis, basis_dr_Y;
        // reverse maps
        Vec< map< Vec<int64_t>, int64_t, vi64less> > cokernels_J_basis_dict, cokernels_I_basis_dict, basis_dr_Y_dict;

        // computes solve_*, cokernels_J_basis*, basis_dr_Y* 
        void init_solve_and_cokernels();
        //  cokernels_I_basis*
        void cokernels_I_basis();

        



        // rho_i,w : P_i --> P_{i-1}
        // m w g \omega/f^{m +1} = w \rho_{i, w} (g) \omega/f^m + m w \pi_i(g) \omega/f^{m+1}
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
        void rho_and_pi_matrices(int64_t d);

        // inclusion_matrices[u] is the matrix that represents:
        // I_u : J_0 + ... + J_n --> P_1 + P_2 + ... + P_{n+1}
        // where ai -- >  u * a_i \in P_{i + 1}
        map< Vec<int64_t>, Mat<int64_t>, vi64less> inclusion_matrices;

        // last_reduction represents
        // M: P_1 + .... + P_{n + 1} --- > P_1 + J_2 + ... + J_n
        // (i - 1)! ai \in P^{i} ---> P_1 + J_2 +... + J_i
        Mat<R> last_reduction;
        //computes the last_reduction
        void last_reduction();


        // projection of matrices 
        // proj_X :  PH^{n-1} (Y) ---> PH^{n-1} (X)
        // proj_notX :  PH^{n-1} (Y) ---> PH^{n-1} (X)^{perp}
        Mat<int64_t> proj_X, proj_notX;
        //computes  proj_X and proj_notX
        void proj();



        
        /*
         * Functions
         */
        
        /*
         * Constructors
         */
        dr(){};
        virtual ~dr(){};
        void init();

        // computes matrix of the map
        // (H0, \dots Hn) ---> h0 * f + \sum_\lambda Hi * \partial_\lambda f
        // where Hi \in P_(d - 1), \lambda runs over a dual basis
        void matrix_J(Mat<R> &result, int64_t d);


        // see solve_system.h
        void solve_system(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<T> &T, const Vec<T> &initB);

        


        
    
};


template<class R>
void dr<R>::compute_f_power(int64_t N)
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
void cokernels_I_basis()
{
    for(int64_t d = 1; d < n + 1; d++)
    {
        Mat<T> MJ;
        matrix_J(MJ, d);
       //FIXME
        //solve system over ZZ_p is broken

    }

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
                result[tuple_dict[d][fit->first + v]][ i * tuple_list[d - 1].length() + j] = fit->second * (dual_basis[i] * fit->first);
    }
}



template<class R>
void dr<R>::solve_system(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<T> &T, const Vec<T> &initB)
{
    if( characteristic != 0)
        solve_system_padic<R>(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<T> &T, const Vec<T> &initB);
    else
        solve_system<R>(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<T> &T, const Vec<T> &initB);
}

#endif
