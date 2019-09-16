// Copyright 2017 Edgar Costa
// See LICENSE file for license details.


#ifndef DR_H
#define DR_H

#include "linear_algebra.h"
#include "tools.h"
#include "finitediff.h"
#include "timing.h"


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

//see void reduce_vector(Vec<R> &H, R &D, const Vec<int64_t> &u, const Vec<int64_t> &v, const int64_t &k, const Vec<R> &G, int64_t method = DEFAULT_VECTOR_REDUCTION_METHOD)
// method = 0, plain
// method = 1, finite diff
// method = 2, finite diff working over ZZ or ZZ[x] to avoid reductions
// method = 3, BSGS not implemented
#define DEFAULT_VECTOR_REDUCTION_METHOD 1

//see void get_reduction_matrix(Vec< Mat<R> > &M, R &Mden, const Vec<int64_t> &u, const Vec<int64_t> &w, int64_t method = DEFAULT_MATRIX_REDUCTION_METHOD);
// method = 0, plain, no pre computations involved
// method = 1, poly, if not already computed, first computes the matrix as polynomial in the entries of  multivariable polynomial
//
// Unless, the dimension and the volume of the polytope are small, I highly recommend to precompute the reduction matrices as multivariable polynomials. There might be cases where this precomputation might become a bottle neck.
#define DEFAULT_MATRIX_REDUCTION_METHOD 1



//R = ZZ, ZZ_p, zz_p. ZZ_pE, or zz_pe
template<typename R>
class dr{
    public:
        /*
         * NOTATION:
         *
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
         * H^{n} (T \ Y) = P_1 + J_2 + ... + J_n = n! vol(Delta) + n
         * PH^{n-1} (Y) = J_1 + J_2 + ... + J_n = n! vol(Delta) - 1
         * PH^{n-1} (X) = P^*_1 + I_2 + ... + I_n
         *
         *
         * NOTE Feb 5, 2018: I think it is H^n(P\X) instead of PH^{n-1}(Y)
         */

        /*
         * Attributes and corresponding initializing functions
         */
        // f_frob = sigma(f)
        map< Vec<int64_t>, R, vi64less> f, f_frob;
        int64_t n;

        /*
         * verbose = 0 --> silent
         * verbose = 1 --> minimal
         * verbose = 2 --> human readable
         * verbose > 2 --> everything
         */
        int64_t verbose;
        bool minimal;

        /*
         * Working ring
         */
        // Specifications that MUST match with the modulus that defines R
        int64_t p; // the characteristic where f comes from, = 0 if R is ZZ, else the characteristic of FFq
        int64_t precision;
        ZZX fE; // 1 if p == q else FFq.defining_polynomial()
        ZZ modulus; // p^precision


        /*
         * Polytope
         *
         * Half space representation of P
         * w \in k*P <=>  AP * v + k *bP >= 0
         */
        Mat<int64_t> AP;
        Vec<int64_t> bP;
        // returns the minimal k such that v \in k*P
        int64_t min_P(const Vec<int64_t> &v){ return ::min_P(AP, bP, v); }
        //returns the minimal k such that v \in k*int(P)
        int64_t min_intP(const Vec<int64_t> &v){ return ::min_intP(AP, bP, v); }



        /*
         * integral points in the polytope
         *
         * tuple_list[d] = integral points in d*P
         * tuple_int_list[d] = integral interior points in d*P
         * we store them as (n + 1) tuples (d, w) where w \in d*P or equivalently w/d \in P
         * d \leq n + 1
         */
        Vec< Vec< Vec<int64_t> > > tuple_list;
        Vec< Vec< Vec<int64_t> > > tuple_int_list;

        // reverse maps
        Vec< map< Vec<int64_t>, int64_t, vi64less> > tuple_dict;
        Vec< map< Vec<int64_t>, int64_t, vi64less> > tuple_int_dict;

        // computes all the above
        void init_tuples();


        // Powers of f
        // f_power[i] = f^i
        // we must use append to extend f_power
        Vec< map< Vec<int64_t>, R, vi64less> > f_power, f_frob_power;
        // if needed computes f^N and adds it to the f_power list
        void init_f_power(int64_t N);
        void init_f_frob_power(int64_t N);



        /*
         *
         * solve_matrix[d] represents the map
         * P_d ---> P_{d-1}^n + J_d
         * g ---> (gi, ci)
         * such that
         * g = \sum gi * Fi + cokernels
         */
        Vec< Mat<R> > solve_matrix;
        Vec< R > solve_denom;

        /*
         * Recall:
         * J_d = (S / J(f))_d
         * I_d = (S^* + J(f) / J(f))_d
         * PH^{n-1} (Y) = P_1 + J_2 + ... + J_n
         * and
         * PH^{n-1} (X) = P^*_1 + I_2 + ... + I_n
         */
        // coKernels_*_dimensions[i] = \dim *_i
        Vec<int64_t> cokernels_J_dimensions, cokernels_I_dimensions;
        int64_t dim_J, dim_I; //sum(coKernels_*_dimensions)

        // cokernels_*_basis[i] = basis for *_i in P_i
        // constructed such that basis for I_i < basis for J_i
        // basis_dr_Y[i] = J_i or P_1 if i == 1
        // basis_dr_X[i] = I_i (note P^*_1 = I_1)
        Vec< Vec< Vec<int64_t> > >cokernels_J_basis, cokernels_I_basis, basis_dr_Y, basis_dr_X;
        // reverse maps
        Vec< map< Vec<int64_t>, int64_t, vi64less> > cokernels_J_basis_dict, cokernels_I_basis_dict, basis_dr_Y_dict, basis_dr_X_dict;
        // \dim PH^{n-1} (Y) and \dim PH^{n-1} (X)
        int64_t dim_dr_Y, dim_dr_X;
        int64_t max_pole; // highest i such that dim J_i != 0 , this happens to be n, unless H^(0,n-1) = 0
        Vec<int64_t> hodge_numbers;


        //  computes cokernels_I_basis* and hodge_numbers
        void init_cokernels_I_basis();
        void init_cokernels_I_basis(int64_t d);

        // computes solve_*, cokernels_J_basis*, basis_dr_*
        // must be called after  init_cokernels_I_basis()
        void init_solve_and_cokernels();





        /*
         * Write
         * rho_i,w : P_i --> P_{i-1}
         * m w g \omega/f^{m +1} = w \rho_{i, w} (g) \omega/f^m + m w \pi_i(g) \omega/f^{m+1}
         * in PH^{n-1}(Y)
         *
         *
         * rho[i] represents rho_i as one degree polynomial in (W_1, ..., W_n) with matrix coefficients
         * rho[i][j] represents the coefficient of rho_i associated to W_{j+1}
         * e.g. rho[i][0] is the constant coefficient
        */
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
        dr(const int64_t &p, const map< Vec<int64_t>, R, vi64less> &f, const map< Vec<int64_t>, R, vi64less> &ffrob, const Mat<int64_t> &AP, const Vec<int64_t> &bP,  const int64_t &verbose = 0, const bool &minimal = false)
        {
            init(p, f, ffrob, AP, bP, verbose, minimal);
        }

        dr(const int64_t &p, const map< Vec<int64_t>, R, vi64less> &f, const Mat<int64_t> &AP, const Vec<int64_t> &bP,  const int64_t &verbose = 0, const bool &minimal = false)
        {
            init(p, f, f, AP, bP, verbose, minimal);
        }
        ~dr(){};
        void init(const int64_t &p, const map< Vec<int64_t>, R, vi64less> &f, const map< Vec<int64_t>, R, vi64less> &ffrob, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const  int64_t &verbose = 0, const bool &minimal = false);

        // computes matrix of the map
        // (H0, \dots Hn) ---> h0 * f + \sum_i Hi * \partial_i f
        // where Hi \in P_(d - 1),  and \partial_i = x_i \partial / \partial x_i
        void matrix_J(Mat<R> &result, int64_t d);

        /*
         * Linear algebra, see linear_algebra.h for more details
         */
        // calls the appropriate version of solve_system(_padic) depending on R
        void solve_system(Vec<int64_t> &B, Mat<R> &Unom, R &Udenom, const Mat<R> &T, const Vec<int64_t> &initB);
        // calls the appropriate version of cokernel_intersection(_local) depending on R
        void cokernel_intersection(Vec<int64_t>  &res, Mat<R> &T, Mat<R> &S);

        /*
         * v \in P1
         * write J_f as J0 + J1 + ... + Jn
         * returns matrix M in R[T] such that represents
         * M: J_0 + J_1 + ... + J_n -->  J_0 + J_1 + ... + J_n
         *    (a0, a1, ..., an) --> (b0, b1, ..., bn)
         * such that
         * u * v^(j + 1) \sum_i (m + j + i)! a_i /f^{m + j + i + 1} \omega
         * = u * v^j \sum_i (m + j + i - 1)! b_i /f^{m + j + i} \omega
         * we also have
         * J_0 --> J_0
         * J_1 --> J_0 + J_1
         * ...
         * J_n --> J_0 + J_1 + ... + J_n
         * where w \in P_m, st (m, j) != 0
         */


        /*
         * computes M
         * by first evaluating RHO at u + T*v
         * and then following the scheme above
         */
        void get_reduction_matrix_plain(Vec< Mat<R> > &M, R &Mden, const Vec<int64_t> &u, const Vec<int64_t> &w);




        /*
         * same as above, but computes the matrix as polynomial in the u coordinates
         * computes a matrix M with coefficients in U0, ..., UN
         * such that M evaluated at the vector u, matches the matrix above
         */
        void compute_reduction_matrix_poly( const Vec<int64_t> &w);
        void get_reduction_matrix_poly(Vec< Mat<R> > &M, R &Mden, const Vec<int64_t> &u, const Vec<int64_t> &w);

        /*
         * where the matrices computed above are stored
         */
        map< Vec<int64_t>, pair<R, map< Vec<int64_t>, Mat<R>, vi64less> >, vi64less> reduction_poly_dict;

        void get_reduction_matrix(Vec< Mat<R> > &M, R &Mden, const Vec<int64_t> &u, const Vec<int64_t> &w, int64_t method = DEFAULT_MATRIX_REDUCTION_METHOD)
        {
            switch(method)
            {
                case 0:
                    {
                        get_reduction_matrix_plain(M, Mden, u, w);
                        break;
                    }
                case 1:
                    {
                        get_reduction_matrix_poly(M, Mden, u, w);
                        break;
                    }
                default:
                    {
                        get_reduction_matrix(M, Mden, u, w);
                        break;
                    }
            }
        }


        /*
         * reduction
         *
         * input:
         * * u \in P_m, m > 0
         * * v \in P_1
         * * k \in NN
         * * G = (a0, a1, ..., an) \in J_0 + \dots + J_n
         * output:
         * * H = (b0, b1, ..., bn) \in J_0 + \dots + J_n
         * * D \in self.R
         * such that
         * D u v^k \sum_i (m + i + k - 1)! ai omega / f^{m + i + k }
         * =
         * u \sum_i (m + i - 1)! bi omega / f^{m + i}
         */
        void reduce_vector(Vec<R> &H, R &D, const Vec<int64_t> &u, const Vec<int64_t> &v, const int64_t &k, const Vec<R> &G, int64_t method = DEFAULT_VECTOR_REDUCTION_METHOD)
        {
            switch(method)
            {
                case 0:
                {
                    reduce_vector_plain(H, D, u, v, k, G);
                    break;
                }
                case 1:
                {
                    reduce_vector_finitediff_plain(H, D, u, v, k, G);
                    break;
                }
                case 2:
                {
                    reduce_vector_finitediff_lift(H, D, u, v, k, G);
                    break;
                }
                default:
                {
                    reduce_vector(H, D, u, v, k, G);
                    break;
                }
            }
        }
        void reduce_vector_plain(Vec<R> &H, R &D, const Vec<int64_t> &u, const Vec<int64_t> &v, const int64_t &k, const Vec<R> &G);
        void reduce_vector_finitediff_plain(Vec<R> &H, R &D, const Vec<int64_t> &u, const Vec<int64_t> &v, const int64_t &k, const Vec<R> &G);
        void reduce_vector_finitediff_lift(Vec<R> &H, R &D, const Vec<int64_t> &u, const Vec<int64_t> &v, const int64_t &k, const Vec<R> &G);
        //TODO
        void reduce_vector_BSGS(Vec<R> &H, const Vec<int64_t> &u, const Vec<int64_t> &v, const int64_t &k, const Vec<R> &G);


        // in a naive way, only using get_reduction_matrix
        // computes the coordinates of x^w / f^m in PH^{n-1} (Y)
        // random = if takes a random path or not
        void monomial_to_basis(Vec<R> &res, R &den, const Vec<int64_t> &w, bool random = false);

        // computes the approximation of Frob(w \omega / f^m) using N terms in the Frobenius expansion
        // Frob(w \omega / f^m) ~ Frob(w) \sum_{j = 0} ^{N - 1} D_{j, m} Frob(f^j) f^{-p(m + j)}

        void frob_monomial(Vec<R> &res, const Vec<int64_t> &w, int64_t N, int64_t method = DEFAULT_VECTOR_REDUCTION_METHOD);

        void frob_matrix(Mat<R> &res, Vec<int64_t> N, int64_t method = DEFAULT_VECTOR_REDUCTION_METHOD);



        /*
         * Test functions
         */


        /*
         * loops over u \in P_k and v \in P_1
         * and checks if the different methods to compute the reduction matrices agree
         */
        void test_reduction_matrices(const int64_t &k);


        /*
         * loops over u \in P_1, and checks if each entry of I_u
         * I_u : J_0 + ... + J_n --> P_1 + P_2 + ... + P_{n+1}
         * corresponds to the inclusion map
         * J_i --> P_{i+1}, where ai -- >  u * a_i \in P_{i + 1}
         */
        void test_inclusion_matrices();
        /*
         * Recall that the last_reduction matrix
         * M: P_1 + .... + P_{n + 1} --- > P_1 + J_2 + ... + J_n
         * wheree  (i - 1)! ai \in P^{i} ---> P_1 + J_2 +... + J_i
         *
         * test_last_reduction()
         * loops over a basis of {u} of PH^{n-1}(Y)
         * and checks if the last_reduction matrix maps
         * u f^k / f^(\deg u + k) to u
         * for all k >= 0 such that \deg u + k < n + 1
         */
        void test_last_reduction();

        /*
         * tests if f^k v \equiv v for v a basis element and k <= N
         */
        void test_monomial_to_basis(int64_t N = 1, bool random = false);

        /*
         * runs all the tests above
         */
        void test_all()
        {
            if( verbose > 0)
                cout << "dr::test_all()" << endl;
            test_last_reduction();
            test_inclusion_matrices();
            test_monomial_to_basis();
            for(int64_t k = 1; k < n + 1; k++)
                test_reduction_matrices(k);

            if( verbose > 0)
                cout << "dr::test_all() done" << endl;
        }

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
    init(local_p, local_f, local_f, local_AP, local_bP, verbose, minimal);
}


template<typename R>
void dr<R>::init(const int64_t &p, const map< Vec<int64_t>, R, vi64less> &f, const map< Vec<int64_t>, R, vi64less> &ffrob, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const int64_t &verbose, const bool &minimal)
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
    v[0] = 1;
    for(fit = ffrob.begin(); fit != ffrob.end(); fit++)
    {
        for(int64_t i = 0; i < n; i++)
            v[i + 1] = fit->first[i];
        this->f_frob[v] = fit->second;
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
        cout<<" verbose = " << this->verbose;
        cout<<endl;
        cout<<"AP = \n" << this->AP<<endl;
        cout<<"bP = " << this->bP<<endl;
        cout << "f = \n";
        cout <<= this->f;
        cout <<endl;
    }

    init_tuples();
    init_cokernels_I_basis();
    f_power.SetLength(0);
    f_frob_power.SetLength(0);
    if(not minimal)
    {
        init_solve_and_cokernels();
        init_rho_and_pi_matrices();
        init_inclusion_matrices();
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
    //changing the length afterwards is troublesome, as maps are not "relocatable"
    tuple_dict.FixAtCurrentLength();
    tuple_int_dict.FixAtCurrentLength();


    //returns false if v \in k*P  and w \notin k*P for some k, i.e., min_P(v) < min_P(w)
    auto minPcompare = [=](const Vec<int64_t> &v, const Vec<int64_t> &w){return ::min_P(AP, bP, v) > ::min_P(AP, bP, w);};



    for(i = 0 ; i < n + 2; i++)
    {
        //we want the integral points to sorted by minPless and then by vi64less
        sort(tuple_int_list[i].begin(), tuple_int_list[i].end(), vi64less() );
        stable_sort(tuple_int_list[i].begin(), tuple_int_list[i].end(), minPcompare);
        stable_sort(tuple_list[i].begin(), tuple_list[i].end(), minPcompare);
        reverse_dict(tuple_dict[i], tuple_list[i]);
        reverse_dict(tuple_int_dict[i], tuple_int_list[i]);
    }

    if(verbose > 1)
    {
        Vec<int64_t> dimP, dimPint;
        dimP.SetLength( n + 2 );
        dimPint.SetLength( n + 2 );
        for(i = 0; i < n + 2; i++)
        {
            dimP[i] = tuple_list[i].length();
            dimPint[i] = tuple_int_list[i].length();
        }
        print(dimP);
        print(dimPint);
    }
    if(verbose > 2)
        cout<<"dr::init_tuples() end"<<endl;
}

template<typename R>
void init_f_power_core(int64_t N, map< Vec<int64_t>, R, vi64less>  &F, Vec< map< Vec<int64_t>, R, vi64less> > &Fpowers )
{
    if(Fpowers.length() < N + 1)
    {
        // maps are not "reallocatable", this avoids the issue
        Vec< map< Vec<int64_t>, R, vi64less> > Fpow;
        //copy Fpowers to Fpow
        Fpow.SetMaxLength(N + 1);
        Fpow.SetLength(Fpowers.length());
        for(int64_t i = 0; i < Fpowers.length(); i++)
            Fpow[i].insert(Fpowers[i].begin(), Fpowers[i].end());

        //deal with the corner case that we have not done anything yet
        if( Fpow.length() < 2)
        {
            Fpow.SetLength(2);
            Vec<int64_t> zero;
            zero.SetLength(F.begin()->first.length(), int64_t(0));
            Fpow[0][zero] = R(1);
            Fpow[1]  = F;
        }
        while(Fpow.length() < N + 1)
        {
            int64_t k = Fpow.length();
            typename map< Vec<int64_t>, R, vi64less>::const_iterator fit, git;
            typename map< Vec<int64_t>, R, vi64less>::iterator itfk;
            map< Vec<int64_t>, R, vi64less> Y;
            Vec<int64_t> u;
            for(git = Fpow[k - 1].cbegin(); git != Fpow[k - 1].cend(); git++)
            {
                for(fit = F.cbegin(); fit != F.cend(); fit++)
                {
                    u = fit->first + git->first;
                    itfk = Y.find(u);
                    if(itfk == Y.end())
                        Y[u] = fit->second * git->second;
                    else
                        itfk->second += fit->second * git->second;
                }
            }
            Fpow.append(Y);
        }
        //replace Fpowers with Fpow (which wasn't reallocated)
        swap(Fpowers, Fpow);
    }
    assert_print(Fpowers.length(), >=, N + 1);
}



template<typename R>
void dr<R>::init_f_power(int64_t N)
{
    if(verbose > 2)
        cout<<"dr::init_f_power("<<N<<")"<<endl;
    init_f_power_core(N, f, f_power);
    if(verbose > 2)
        cout<<"dr::init_f_power("<<N<<") done"<<endl;
}
template<typename R>
void dr<R>::init_f_frob_power(int64_t N)
{
    if(verbose > 2)
        cout<<"dr::init_f_power("<<N<<")"<<endl;
    init_f_power_core(N, f_frob, f_frob_power);
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
    cokernels_J_dimensions.SetLength(n + 2, int64_t(0));
    assert_print(coeff(P, n + 1), ==, 0);
    for(i = 0; i < n + 2; i++)
    {
        cokernels_J_dimensions[i] = conv< long >( coeff(P, i) );
        if(cokernels_J_dimensions[i] == 0)
        {
            this->max_pole = i - 1;
            break;
        }
    }


    dim_J = sum(cokernels_J_dimensions);

    if( verbose > 1 )
        cout << "dim J_i = "<<cokernels_J_dimensions<<endl;

    cokernels_J_basis.SetLength(max_pole + 2);
    cokernels_J_basis_dict.SetLength(max_pole + 2);
    //changing the length afterwards is troublesome, as maps are not "relocatable"
    cokernels_J_basis_dict.FixAtCurrentLength();

    solve_matrix.SetLength(max_pole + 2);
    solve_denom.SetLength(max_pole + 2);
    for(i =0; i <= max_pole + 1; i++)
    {
        if(i == max_pole + 1 and verbose > 1)
            cout << "Asserting that the hypersurface is non degenerate" << endl;

        Mat<R> J;
        matrix_J(J, i);
        if( verbose > 1)
            printf("Solving Jacobian relations at degree %ld (%ld x %ld)\n", (long)i, (long)J.NumRows(), (long)J.NumCols());

        Vec<int64_t> B, initB;
        if( i <= max_pole )
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
            cout<<"Something went wrong, perhaps the hypersurface is not nondegenerate..."<<endl;
            printf("At degree %ld expected B.length() = %ld, got %ld", (long) i,  (long) cokernels_J_dimensions[i], (long) B.length());
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

    // PH^{n-1} (Y) = P_1 + J_2 + ... + J_{max_pole}, generally max_pole = n
    // PH^{n-1} (X) = P^*_1 + I_2 + ... + I_{max_pole}
    basis_dr_Y.SetLength(max_pole + 1);
    basis_dr_Y_dict.SetLength(max_pole + 1);


    basis_dr_X.SetLength(max_pole + 1);
    basis_dr_X_dict.SetLength(max_pole + 1);

    //changing the length afterwards is troublesome, as maps are not "relocatable"
    basis_dr_Y_dict.FixAtCurrentLength();
    basis_dr_X_dict.FixAtCurrentLength();


    basis_dr_Y[0].SetLength(0);
    basis_dr_Y[1] = tuple_list[1];
    basis_dr_Y_dict[1] = tuple_dict[1];
    dim_dr_Y = tuple_list[1].length();

    basis_dr_X[0].SetLength(0);
    basis_dr_X[1] = tuple_int_list[1];
    basis_dr_X_dict[1] = tuple_int_dict[1];
    dim_dr_X =  tuple_int_list[1].length();


    for(i = 2; i < max_pole + 1; i++)
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
    //FIXME ?? I don't remember what I was supposed to fix
    if(verbose > 1)
        cout<<"dr::init_cokernels_I_basis()"<<endl;
    cokernels_I_dimensions.SetLength(n + 1, int64_t(0));
    cokernels_I_basis.SetLength(n + 1);
    cokernels_I_basis_dict.SetLength(n + 1);
    //changing the length afterwards is troublesome, as maps are not "relocatable"
    cokernels_I_basis_dict.FixAtCurrentLength();
    hodge_numbers.SetLength(n, int64_t(0));
    hodge_numbers[0] = tuple_int_list[1].length();
    int64_t max_hi = n;
    for(int64_t d = 1; d < max_hi + 1; d++)
    {
        init_cokernels_I_basis(d);
        reverse_dict(cokernels_I_basis_dict[d], cokernels_I_basis[d]);
        cokernels_I_dimensions[d] = cokernels_I_basis[d].length();
        hodge_numbers[d - 1] = cokernels_I_basis[d].length();
        if( hodge_numbers[d - 1] == 0 and 2*d <  n + 1 )
        {
            // hodge symmetry implies H^(w -d, d) = H^(d, w-d) = 0
            max_hi--;
            hodge_numbers[n - 1 - (d - 1)] = 0;
            cokernels_I_basis[n - d + 1].SetLength(0);
        }
        if( 2*d >  n + 1) // d - 1 > (n - 1)- ( d - 1)
        {
            if(hodge_numbers[d - 1] != hodge_numbers[(n - 1) - (d -1)] )
            {
                cout<<"Something went wrong, perhaps the hypersurface is not nondegenerate..."<<endl;
                print(d);
                print(n - 1 - (d-1));
                assert_print(hodge_numbers[d - 1], ==, hodge_numbers[(n - 1) - (d -1)]);
            }
        }
    }
    dim_I = sum(cokernels_I_dimensions);
    if(verbose > 1)
    {
        cout << "dim I_i = "<<cokernels_I_dimensions<<endl;
        print(hodge_numbers);
    }
    if(verbose > 2)
        cout<<"dr::init_cokernels_I_basis() done"<<endl;


}

template<typename R>
void dr<R>::init_cokernels_I_basis(int64_t d)
{
    if(verbose > 2)
        cout<<"dr::init_cokernels_I_basis("<<d<<")"<<endl;
    if(d == 1)
    {
        cokernels_I_basis[d] = tuple_int_list[d];
    }
    else
    {
        //I_d = (S^* + J(f) ) / J(f))_d
        Mat<R> J;
        Vec<int64_t> nonpivots;
        int64_t i, dim_Sintd;
        matrix_J(J, d);
        dim_Sintd = tuple_int_list[d].length();

        //Matrix representing S^*_d in S_d
        Mat<R> T;
        T.SetDims(J.NumRows(), dim_Sintd);
        for(i = 0; i < dim_Sintd; i++)
            T[tuple_dict[d][tuple_int_list[d][i]]][i] = 1;

        //computes the cokernel of the map (img(J) \cap S^* )_d -> S^*_d
        cokernel_intersection(nonpivots, T, J);
        cokernels_I_basis[d].SetLength(nonpivots.length());
        for(i = 0; i < nonpivots.length(); i++)
            cokernels_I_basis[d][i] = tuple_int_list[d][nonpivots[i]];
    }
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
    timestamp_type time1, time2;
    get_timestamp(&time1);
    ::solve_system_padic<R>(B, Unom, Udenom, T, initB, p, precision, fE);
    get_timestamp(&time2);
    if( verbose > 1)
        printf("Time elapsed solve_system_padic (%lu x %lu matrix) %f s\n",(long unsigned)T.NumRows(), (long unsigned)T.NumCols(), timestamp_diff_in_seconds(time1,time2));
}

template <> inline void dr<ZZ>::solve_system(Vec<int64_t> &B, Mat<ZZ> &Unom, ZZ &Udenom, const Mat<ZZ> &T, const Vec<int64_t> &initB)
{
    timestamp_type time1, time2;
    get_timestamp(&time1);
    ::solve_system(B, Unom, Udenom, T, initB);
    get_timestamp(&time2);
    if( verbose > 1)
        printf("Time elapsed solve_system (%lu x %lu matrix) over ZZ %f s\n",(long unsigned)T.NumRows(), (long unsigned)T.NumCols(), timestamp_diff_in_seconds(time1,time2));
}

template<typename R>
void dr<R>::cokernel_intersection(Vec<int64_t>  &res, Mat<R> &T, Mat<R> &S)
{
    timestamp_type time1, time2;
    get_timestamp(&time1);
    ::cokernel_intersection_local(res, T, S, p, fE);
    get_timestamp(&time2);
    if( verbose > 1)
        printf("Time elapsed cokernel_intersection_local (%lu x %lu matrix) %f s\n",(long unsigned)T.NumRows(), (long unsigned)T.NumCols(), timestamp_diff_in_seconds(time1,time2));
    }
template<>
inline void dr<ZZ>::cokernel_intersection(Vec<int64_t>  &res, Mat<ZZ> &T, Mat<ZZ> &S)
{
    timestamp_type time1, time2;
    get_timestamp(&time1);
    ::cokernel_intersection(res, T, S);
    get_timestamp(&time2);
    if( verbose > 1)
        printf("Time elapsed cokernel_intersection (%lu x %lu matrix) over ZZ %f s\n",(long unsigned)T.NumRows(), (long unsigned)T.NumCols(), timestamp_diff_in_seconds(time1,time2));
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
    rho.SetLength(max_pole + 2);
    rho_den.SetLength(max_pole + 2);
    pi.SetLength(max_pole + 2);
    pi_den.SetLength(max_pole + 2);

    //deal with d = 0
    // rho[0] = 0 map, the codomain doesn't exist...
    // pi[0] is the identity on {0}
    rho[0].SetLength(0);
    rho_den[0] = 1;
    pi[0].SetDims(1, 1);
    pi[0][0][0] = 1;
    pi_den[0] = 1;


    for(d = 1; d < max_pole + 2; d++)
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
    for(int64_t k = 0; k < max_pole + 1; k++)
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
        for(int64_t k = 0; k < max_pole + 1; k++)
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
    /*
     * loops over u \in P_1, and checks if each entry of I_u
     * I_u : J_0 + ... + J_n --> P_1 + P_2 + ... + P_{n+1}
     * corresponds to the inclusion map
     * J_i --> P_{i+1}, where ai -- >  u * a_i \in P_{i + 1}
     */

    int64_t i, j, k;
    Vec< Vec<int64_t> > B, BJ;
    for(i = 1; i < max_pole + 2; i++)
        B.append(tuple_list[i]);
    for(i = 0; i < max_pole + 1; i++)
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
                    assert_print(M[i][j], ==, 1);
                    assert_print(u + wj, ==, B[i]);
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

    shift_rows.SetLength(n + 2, int64_t(0));
    shift_columns.SetLength(n + 2, int64_t(0));
    // shift_rows[k] = dim(P_1 + J_2 + ... + J_{k - 1})
    // shift_columns[k] = dim(P_1 + P_2 + .... + P_{k - 1}
    int64_t total_columns, total_rows;
    total_columns = 0;
    total_rows = 0;
    for(k = 1; k < max_pole + 2; k++)
    {
        shift_columns[k] = total_columns;
        shift_rows[k] = total_rows;
        total_columns += tuple_list[k].length();
        if(k < max_pole + 1)
            total_rows += basis_dr_Y[k].length();
    }
    M.SetDims(total_rows, total_columns);

    D = 1;
    for(k = 2; k < max_pole + 2; k++)
        D *= pi_den[k];

    //P_1 --> P_1 is the identity map
    for(k = 0; k < tuple_list[1].length(); k++)
        M[k][k] = D;


    int64_t factorial = 1;
    for(k = 2; k < max_pole + 2; k++)
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

    /*
     * Recall that the last_reduction matrix
     * M: P_1 + .... + P_{n + 1} --- > P_1 + J_2 + ... + J_n
     * wheree  (i - 1)! ai \in P^{i} ---> P_1 + J_2 +... + J_i
     *
     * test_last_reduction()
     * loops over a basis of {u} of PH^{n-1}(Y)
     * and checks if the last_reduction matrix maps
     * u f^k / f^(\deg u + k) to u
     * for all k >= 0 such that \deg u + k < n + 1
     */
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
        for(k = 0; k < max_pole + 2 - v[0]; k++)
        {
            Vec<R> G;
            G.SetLength(last_reduction.NumCols(), R(0));
            shift = 0;
            for(j = 1; j < k + v[0]; j++)
                shift += tuple_list[j].length();

            typename map< Vec<int64_t>, R, vi64less>::const_iterator fkit;
            for(fkit = f_power[k].cbegin(); fkit != f_power[k].cend(); fkit++)
                G[shift + tuple_dict[k + v[0]][fkit->first + v]] += fkit->second;

            if(  factorial<R>(k + v[0] - 1) * expected != last_reduction * G)
            {
                cout << "TEST FAILED" <<endl;
                print(i);
                print(k);
                print(v);
                print(G);
                print(factorial<R>(k + v[0] - 1));
                print(expected);
                print(factorial<R>(k + v[0] - 1) * expected);
                print(last_reduction * G);
                abort();
            }
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
    for(i = 1; i < max_pole + 1; i++)
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
void dr<R>::get_reduction_matrix_plain(Vec< Mat<R> > &M, R &Mden, const Vec<int64_t> &u, const Vec<int64_t> &v)
{
    if( verbose > 2)
        cout << "dr::get_reduction_matrix(u = "<<u<<", v = "<<v<<")" << endl;

    //asserts \deg v == 1
    assert_print(v[0], ==, 1);

    int64_t i, j, k, l, z;


    // M = M[0] + M[1] * T + ... + M[n] * T^(n + 1)
    M.kill();
    M.SetLength(max_pole + 2);
    for(i = 0; i < max_pole + 2; i++)
        M[i].SetDims(dim_J, dim_J);

    // rho_{k, u + T*v} = RHO0[k] + T*RHO1[k] 1 <= k <= n + 1
    Vec< Mat<R> > RHO0, RHO1;;
    RHO0.SetLength(max_pole + 2);
    RHO1.SetLength(max_pole + 2);
    Mden = R(1);
    for(k = 1; k < max_pole + 2; k++)
    {
        RHO0[k] = rho[k][0];
        RHO1[k].SetDims(rho[k][0].NumRows(), rho[k][0].NumCols());
        for(i = 0; i < n + 1; i++)
        {
            RHO0[k] += u[i] * rho[k][i + 1];
            RHO1[k] += v[i] * rho[k][i + 1];
        }
        Mden *= rho_den[k];
    }
    Vec<int64_t> shift; //shift[k] = sum(map(len, cokernels_J_basis[:k]))
    shift.SetLength(max_pole + 1, int64_t(0));
    for(k = 1; k < max_pole + 1; k++)
    {
        shift[k] = shift[k-1] + cokernels_J_basis[k-1].length();
    }

    /*
     * let alpha be a basis element in J_k
     * M(alpha) = (c_0, c_1, c_2, ..., c_k, 0, ...)
     * where
     * h_{k + 1} = alpha + v \in P_{k + 1}
     * c_i = pi_{i} (h_i)  \in J_i  with i <= k
     * h_i = rho_{i + 1) u + Y*v} (h_{i + 1}) \in P_i with i <= k
     */
    for(k = 0; k < max_pole + 1; k++)
    {
        for(z = 0; z < cokernels_J_dimensions[k]; z++)
        {
            Vec<int64_t> &b = cokernels_J_basis[k][z];
            int64_t b_coordinate = shift[k] + z; // as an element in cokernels_J_basis

            // hi is represented as a list [hi0, hi1, ... ]
            // where hi = \sum_k hik Y^k
            // h_{k+1} \in P_{k + 1}
            Vec< Vec<R> > hi;
            hi.SetLength(1);
            hi[0].SetLength( tuple_list[k+1].length() );
            hi[0][ tuple_dict[k + 1][b + v] ] = Mden;
            for(i = k + 1; i >= 0; i--)
            {
                assert_print(hi.length(), ==, (k + 1 - i) + 1);

                //copy pi_i(hi) to the matrix
                for(j = 0; j < hi.length(); j++)
                {
                    //hi \in P_i
                    // ci = pi_i(h_i) \in J_i
                    // ci*T^j corresponds to pi_i(h_i*T^j)
                    // i.e., the jth entry in hi and M
                    Vec<R> cij;
                    cij = pi[i] * hi[j];
                    if( pi_den[i] != R(1) )
                        cij = cij/pi_den[i];

                    for(l = 0; l < cokernels_J_dimensions[i]; l++)
                        M[j][shift[i] + l][b_coordinate] += cij[l];
                }

                if(i > 0)
                {
                    //h_{i - 1} = \rho_{i} (hi) \in P_{i - 1}
                    Vec< Vec<R> > hnew;
                    hnew.SetLength( hi.length() + 1 );
                    for(l = 0; l < hnew.length(); l++)
                        hnew[l].SetLength( tuple_list[i - 1].length() );
                    for(l = 0; l < hi.length(); l++)
                    {
                        hnew[l] += RHO0[i] * hi[l];
                        hnew[l + 1] += RHO1[i] * hi[l];
                    }
                    if( rho_den[i] != R(1) )
                        for(l = 0; l < hnew.length(); l++)
                            hnew[l] = hnew[l]/rho_den[i];

                    hi = hnew;
                }

            }
        }
    }


    if( verbose > 2)
        cout << "dr::get_reduction_matrix(u = "<<u<<", v = "<<v<<") done" << endl;
}


template<typename R>
void dr<R>::compute_reduction_matrix_poly(const Vec<int64_t> &v)
{
    // do we really need to compute it?
    typename map< Vec<int64_t>, pair< R, map< Vec<int64_t>, Mat<R>, vi64less> >, vi64less>::const_iterator it;
    it = reduction_poly_dict.find(v);
    if( it != reduction_poly_dict.end())
        return;

    if( verbose > 1)
        cout << "dr::compute_reduction_poly(v = "<<v<<")\n";
    timestamp_type time1, time2;
    double wall_time, user_time;
    user_time = get_cpu_time();
    get_timestamp(&time1);

    //asserts \deg v == 1
    assert_print(v[0], ==, 1);

    int64_t i, j, k, l , z;


    map< Vec<int64_t>, Mat<R>, vi64less> &M = reduction_poly_dict[v].second;
    R &Mden = reduction_poly_dict[v].first;
    Mden = R(1);



    Vec<int64_t> shift; //shift[k] = sum(map(len, cokernels_J_basis[:k]))
    shift.SetLength(max_pole + 2, int64_t(0));
    for(k = 1; k < max_pole + 2; k++)
        shift[k] = shift[k-1] + cokernels_J_basis[k-1].length();

    for(k = 1; k < max_pole + 2; k++)
        Mden *= rho_den[k];


    Vec<int64_t> zero;
    zero.SetLength(n + 1, 0);


    // let alpha be a basis element in J_k
    // M(alpha) = (c_0, c_1, c_2, ..., c_k, 0, ...)
    // where
    // h_{k + 1} = alpha + v \in P_{k + 1}
    // c_i = pi_{i} (h_i)  \in J_i  with i <= k
    // h_i = rho_{i + 1) u + Y*v} (h_{i + 1}) \in P_i with i <= k
    for(k = 0; k < max_pole + 2; k++)
    {
        for(z = 0; z < cokernels_J_dimensions[k]; z++)
        {
            Vec<int64_t> &b = cokernels_J_basis[k][z];
            int64_t b_coordinate = shift[k] + z; // as an element in cokernels_J_basis

            // hi \in Vec(R[W0, ..., WN]) is represented as a dictionary
            // monomial --> coefficient
            map< Vec<int64_t>, Vec<R>, vi64less> hi;
            hi[zero].SetLength( tuple_list[k + 1].length() );
            hi[zero][ tuple_dict[k + 1][b + v] ] = Mden;
            for(i = k + 1; i >= 0 ; i-- )
            {
                //deg of monomial is k + 1 - i, as expected
                assert_print( (hi.begin()->first)[0], <= , k + 1 - i);

                // copy pi_i(hi) to the matrix
                typename map< Vec<int64_t>, Vec<R>, vi64less>::const_iterator cit;
                for(cit = hi.begin(); cit != hi.end(); cit++)
                {
                    // ci = pi_i(hi) \in J_i
                    // ci[v] corresponds to pi_i( h_i[v])
                    // where u = zero or ej
                    const Vec<int64_t> &w = cit->first;
                    const Vec<R> &hiw = cit->second;
                    Vec<R> ciw;
                    ciw = pi[i] * hiw;
                    if( pi_den[i] != R(1) )
                        ciw = ciw/pi_den[i];

                    Mat<R> &Mw = M[w];
                    Mw.SetDims(dim_J, dim_J);
                    for(l = 0; l < cokernels_J_dimensions[i]; l++)
                    {
                        Mw[shift[i] + l][b_coordinate] += ciw[l];
                    }
                }
                if(i > 0)
                {
                    //h_{i -1} = \rho_{i} (hi) \in P_{i-1}
                    map< Vec<int64_t>, Vec<R>, vi64less> hnew;
                    int64_t hnew_length = tuple_list[i - 1].length();
                    for(cit = hi.begin(); cit != hi.end(); cit++)
                    {
                        Vec<int64_t> w = cit->first;
                        const Vec<R> &hiw = cit->second;

                        // recall rho[i][0] corresponds to the constant coefficient
                        // and rho[i][j+1] corresponds to the W_j coefficient
                        //deal with the constant coefficient
                        if( hnew.find(w) == hnew.end() )
                            hnew[w].SetLength(hnew_length, conv<R>(0));
                        hnew[w] += rho[i][0] * hiw;

                        for(j = 0; j < n + 1; j++)
                        {
                            Vec<int64_t> t = w;
                            t[j]++;
                            // assure the right length
                            if( hnew.find(t) == hnew.end() )
                                hnew[t].SetLength(hnew_length, conv<R>(0));
                            hnew[t] += rho[i][j+1] * hiw;
                        }
                    }

                    // deal with denominators
                    if( rho_den[i] != R(1) )
                    {
                        typename map< Vec<int64_t>, Vec<R>, vi64less>::iterator it;
                        for( it = hnew.begin(); it != hnew.end(); it++ )
                            it->second = it->second / rho_den[i];
                    }
                    hi.swap(hnew);
                }
            }
        }
    }
    get_timestamp(&time2);
    wall_time = timestamp_diff_in_seconds(time1,time2);
    user_time = get_cpu_time() - user_time;
    if( verbose > 1 )
    {
        cout << "dr::compute_reduction_poly(v = "<<v<<") done ";
        printf("Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );
    }
}

template<typename R>
void dr<R>::get_reduction_matrix_poly(Vec< Mat<R> > &M, R &Mden, const Vec<int64_t> &u, const Vec<int64_t> &v)
{
    if( verbose > 2)
        cout << "dr::get_reduction_matrix(u = "<<u<<", v = "<<v<<")" << endl;
    // timestamp_type time1, time2;
    // double wall_time, user_time;
    // user_time = get_cpu_time();
    // get_timestamp(&time1);

    typename map< Vec<int64_t>, pair< R, map< Vec<int64_t>, Mat<R>, vi64less> > , vi64less>::const_iterator it;

    it = reduction_poly_dict.find(v);
    if( it == reduction_poly_dict.end())
    {
        compute_reduction_matrix_poly(v);
        it = reduction_poly_dict.find(v);
        assert(it != reduction_poly_dict.end());
    }
    //Now we evaluate the polynomial matrix M(U0, U1, .., UN) at M(u + Y v)
    Mden = (it->second).first;
    const map< Vec<int64_t>, Mat<R>, vi64less> &Mdict = (it->second).second;

    int64_t i;

    // initialize M
    // release space and set to length 0
    M.kill();
    M.SetLength(max_pole + 2);
    for(i = 0; i < max_pole + 2; i++)
        M[i].SetDims(dim_J, dim_J);

    typename map< Vec<int64_t>, Mat<R>, vi64less>::const_iterator itM;
    ZZX monomial_evaluated;
    Vec<ZZ> monomial;
    for(itM = Mdict.begin(); itM != Mdict.end(); itM++)
    {
        // itM->first = alpha = (alpha0, alpha1,..., alpha_n) representing U^alpha
        // we want to compute (u + Y*v)^alpha
        const Vec<int64_t> &alpha = itM->first;
        //monomial_evaluated = ZZX(1);
        set(monomial_evaluated);
        for( i = 0; i < n + 1; i++)
        {
            if(alpha[i] != 0 and u[i] == 0 and v[i] == 0)
            {
                clear(monomial_evaluated);
                break;
            }
            binomial_expansion(monomial, u[i], v[i], alpha[i]);
            monomial_evaluated *= conv<ZZX>(monomial);
        }
        assert_print( deg(monomial_evaluated), <= , max_pole + 1);
        for(i = 0; i < deg(monomial_evaluated) + 1; i++)
            M[i] += conv<R>(monomial_evaluated[i]) * (itM->second);
    }

    //get_timestamp(&time2);
    //wall_time = timestamp_diff_in_seconds(time1,time2);
    //user_time = get_cpu_time() - user_time;
    if( verbose > 2)
    {
        cout << "dr::get_reduction_matrix(u = "<<u<<", v = "<<v<<") done" << endl;
        //printf("Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );
    }
}


template<typename R>
void dr<R>::test_reduction_matrices(const int64_t &k)
{
    int64_t i, j, l;
    for(l = 0; l < tuple_list[1].length(); l++)
    {
        Vec<int64_t> &v = tuple_list[1][l];
        for(i = 0; i < tuple_list[k].length(); i++)
        {
            Vec<int64_t> &u = tuple_list[k][i];
            Vec< Mat<R> > Mpoly, Mplain;
            R den_poly, den_plain;
            get_reduction_matrix_poly( Mpoly, den_poly, u, v);
            get_reduction_matrix_plain( Mplain, den_plain, u, v);
            for(j = 0; j < max_pole + 1; j++)
            {
                if( Mplain[j] * den_poly != Mpoly[j] * den_plain)
                {
                    cout << "FAIL" <<endl;
                    print(Mplain.length());
                    print(u);
                    print(v);
                    print(j);
                    cout << endl;
                    print(Mplain[j]);
                    print(den_plain);
                    cout << endl;
                    print(Mpoly[j]);
                    print(den_poly);
                    print(Mplain[j] - Mpoly[j]);

                    cout << endl << endl;
                    cout <<= reduction_poly_dict[v].second;
                    cout << endl;

                    assert(false);
                }
            }
        }
    }
}


template<typename R>
void dr<R>::reduce_vector_plain(Vec<R> &H, R &D, const Vec<int64_t> &u, const Vec<int64_t> &v, const int64_t &k, const Vec<R> &G)
{
    if( verbose > 2) {
        cout << "dr::reduce_vector_plain(u = "<<u<<", v = "<<v<<", k = "<<k<<")" << endl;
        cout << G << endl;
    }

    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
    user_time = get_cpu_time();
    get_timestamp(&wtime1);

    int64_t dim_low = dim_J - cokernels_J_dimensions[n]; // sum(self.coKernels_J_dimensions[:self.n])
    Vec< Mat<R> > M;
    R Mden;
    get_reduction_matrix(M, Mden, u, v);
    Vec<R> Gin, Gout;
    Gin = G;
    D = 1;
    double ttime, utime;
    ttime = 0;
    for(int64_t l = k - 1; l > 0; l--)
    {
        D *= Mden;
        R lpower = R(1);
        Gout = M[0] * Gin;
        utime = get_cpu_time();
        for(int64_t i = 1; i < n + 1; i++)
        {
            lpower *= l;//lpower = l^i
            Gout += lpower * (M[i] * Gin);
        }
        ttime += get_cpu_time() - utime;
        //taking advantage that \dim J_0 = 1
        lpower *= l; //lpower = l^(n+1)
        for(int64_t i = 0; i < cokernels_J_dimensions[n]; i++)
            Gout[0] += lpower * (M[n + 1][0][dim_low + i] * Gin[dim_low + i]);

        Gin = Gout;
    }
    H = M[0] * Gin;
    D *= Mden;
    if( verbose > 2)
    {
        get_timestamp(&wtime2);
        wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
        user_time = get_cpu_time() - user_time;
        cout << "dr::reduce_vector_plain(u = "<<u<<", v = "<<v<<") done";
        printf(" Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );
    }
}

template<typename R>
void dr<R>::reduce_vector_finitediff_plain(Vec<R> &H,  R &D, const Vec<int64_t> &u, const Vec<int64_t> &v, const int64_t &k, const Vec<R> &G)
{
    if( verbose > 2)
        cout << "dr::reduce_vector_finitediff_plain(u = "<<u<<", v = "<<v<<", k = "<<k<<")" << endl;
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
    user_time = get_cpu_time();
    get_timestamp(&wtime1);

    if(k <= n +1)
        reduce_vector_plain(H, D, u, v, k, G);
    else
    {
        Vec< Mat<R> > M;
        R Mden;
        get_reduction_matrix(M, Mden, u, v);
        finitediff_plain(H, k, G, M);
        D = power(Mden, k);
    }
    if( verbose > 2)
    {
        get_timestamp(&wtime2);
        wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
        user_time = get_cpu_time() - user_time;
        cout << "dr::reduce_vector_finitediff_plain(u = "<<u<<", v = "<<v<<") done";
        printf(" Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );
    }
}
template<typename R>
void dr<R>::reduce_vector_finitediff_lift(Vec<R> &H,  R &D, const Vec<int64_t> &u, const Vec<int64_t> &v, const int64_t &k, const Vec<R> &G)
{
    if( verbose > 2)
        cout << "dr::reduce_vector_finitediff_lift(u = "<<u<<", v = "<<v << ", k = "<<k<<")" << endl;
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
    user_time = get_cpu_time();
    get_timestamp(&wtime1);

    if(k <= n +1)
        reduce_vector_plain(H, D, u, v, k, G);
    else
    {
        Vec< Mat<R> > M;
        R Mden;
        get_reduction_matrix(M, Mden, u, v);
        finitediff_lift(H, k, G, M);
        D = power(Mden, k);
    }
    if( verbose > 2)
    {
        get_timestamp(&wtime2);
        wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
        user_time = get_cpu_time() - user_time;
        cout << "dr::reduce_vector_finitediff_lift(u = "<<u<<", v = "<<v<<") done";
        printf(" Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );
    }
}



template<typename R>
void dr<R>::monomial_to_basis(Vec<R> &G, R &D, const Vec<int64_t> &w, bool random)
{
    if(verbose > 1)
        cout << "dr::monomial_to_basis(w = "<<w<<")" << endl;
    int64_t e, i, m;
    Vec<int64_t> u, v;
    // w \in P_m
    m = w[0];
    u = w;

    G.kill();
    G.SetLength(dim_J, R(0));
    G[0] = 1;
    D = 1;


    // W = (m-1)! u G / f^m \omega
    for(e = m; e > 1; e--)
    {
        assert_print(e, ==, u[0]);
        if( random )
        {
            while(true)
            {
                v =  tuple_list[1][rand() % tuple_list[1].length()];
                if( min_P(u - v) <= e - 1)
                    break;
            }
        }
        else
        {
            for(i = 0; i < tuple_list[1].length(); i++)
            {
                v = tuple_list[1][i];
                if( min_P(u - v) <= e - 1)
                    break;
            }
        }
        u = u - v;
        if(verbose > 2)
            cout << "\t" << u << " + " << v <<" --> "<<u<<endl;
        Vec< Mat<R> > M;
        R Mden;
        get_reduction_matrix(M, Mden, u, v);
        G = M[0]*G;
        D *= Mden;
    }
    assert_print(u[0], ==, 1);
    G = last_reduction * (inclusion_matrices[tuple_dict[1][u]] * G);
    D *= last_reduction_den;

    if( m > 0 )
        D *= factorial<R>(m - 1);
    if(verbose > 1)
        cout << "dr::monomial_to_basis(w = "<<w<<") done" << endl;
}


template<typename R>
void dr<R>::test_monomial_to_basis(int64_t N, bool random)
{
    if( verbose > 1 )
        cout << "dr::test_monomial_to_basis()" << endl;

    Vec< Vec<int64_t> > B;
    int64_t i, k;

    for(i = 0; i < basis_dr_Y.length(); i++)
        B.append( basis_dr_Y[i] );

    init_f_power(N);
    for(k = 0; k < N + 1; k++)
    {
        if(verbose > 2)
            cout << "testing if  v f^"<<k<<" == v"<<endl;

        for(i = 0; i < B.length(); i++)
        {
            Vec<int64_t> &v = B[i];
            Vec<R> expected, res;
            expected.SetLength(B.length(), R(0));
            res.SetLength(B.length(), R(0));
            expected[i] = 1;
            typename map< Vec<int64_t> , R, vi64less>::const_iterator fit;
            R resden;
            resden = 0;
            for(fit = f_power[k].cbegin(); fit != f_power[k].cend(); fit++)
            {
                Vec<R> G;
                R D;
                monomial_to_basis(G, D, fit->first + v, random);
                if(resden == 0)
                    resden = D;
                else
                    assert_print(resden, ==, D);

                res += fit->second * G;
            }

            if( resden * expected !=  res)
            {
                cout << "TEST FAILED" <<endl;
                print(k);
                print(v);
                print(expected);
                print(res);
                print(resden);
                abort();
            }
        }
    }
    if( verbose > 2 )
        cout << "dr::test_monomial_to_basis() done" << endl;
}


template<typename R>
void dr<R>::frob_matrix(Mat<R> &res, Vec<int64_t> N, int64_t method)
{
    if(verbose > 0)
        cout <<"dr<R>::frob_matrix("<<N<<") "<<endl;
    assert_print( max_pole, <=, N.length() );
    Mat<R> F;
    int64_t i, m, shift;
    F.kill();
    F.SetDims( dim_dr_X, dim_dr_Y );
    shift = 0;
    for(m = 1; m < max_pole + 1; m++)
    {
        for(i = 0; i < basis_dr_X[m].length(); i++)
        {
            if(verbose > 0)
                cout<<"Computing Frob("<<basis_dr_X[m][i]<<") m = "<<m<<" N = "<<N[m-1]<<endl;
            frob_monomial(F[shift + i], basis_dr_X[m][i], N[m - 1], method);
        }
        shift += basis_dr_X[m].length();
    }
    F = transpose(F);
    res = proj_X * F;

    Mat<R> zero;
    zero.SetDims( dim_dr_Y, dim_dr_X);
    assert_print(zero, ==, proj_notX * F);
    if(verbose > 0)
        cout <<"dr<R>::frob_matrix("<<N<<") done!"<<endl;
}

template<typename R>
void dr<R>::frob_monomial(Vec<R> &F, const Vec<int64_t> &w, int64_t N, int64_t method)
{
    if(verbose > 0)
        cout <<"dr<R>::frob_monomial("<<w<<", "<<N<<") "<<endl;
    if(N == 0)
    {
        if(F.length() != dim_dr_Y)
            F.SetLength(dim_dr_Y);
        for(int64_t i = 0; i < dim_dr_Y; i++)
            F[i] = R(0);
        if(verbose > 0)
            cout <<"dr<R>::frob_monomial("<<w<<", "<<N<<") done!"<<endl;
        return;
    }
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
    user_time = get_cpu_time();
    get_timestamp(&wtime1);

    assert_print(p, !=, 0);

    init_f_frob_power(N - 1);
    int64_t i, j, m, e, end;
    R fact = R(1);

    map< Vec<int64_t>, Vec<R>, vi64less> H;
    typename map< Vec<int64_t>, Vec<R>, vi64less>::const_iterator Hit;
    typename map< Vec<int64_t>, Vec<R>, vi64less>::iterator Hnewit;
    m = w[0];

    for(e = m + N - 1; e > 0; e--)
    {
        map< Vec<int64_t>, Vec<R>, vi64less> Hnew;
        if( verbose > 1 )
            cout <<"\t e = "<< e << endl;

        //add terms of pole order e
        if( e >= m )
        {
            //binomial(-m, j) = (-1) ^i binomial(m + i -1, m - 1)
            // Djm = binomial(-m, j) * binomial(m + N - 1, N -j - 1)
            j = e - m;
            R Djm;
            Djm = conv<R>( binomial(m + j - 1, m - 1));
            if (j%2 != 0)
                Djm = -Djm;
            Djm *= conv<R>( binomial(m + N - 1, N -j - 1) );

            typename map< Vec<int64_t>, R, vi64less >::const_iterator fjit;
            for(fjit = f_frob_power[j].cbegin(); fjit != f_frob_power[j].cend(); fjit++)
            {
                if(fjit->second != 0)
                {
                    Vec<int64_t> monomial = w + fjit->first;
                    Hit = H.find(monomial);
                    if( Hit == H.end() )
                    {
                        H[monomial].SetLength(dim_J, R(0));
                    }
                    H[monomial][0] += Djm * fact * fjit->second;
                }
            }
        }

        end = (e > 1) ? p : p - 1;

        for(i = 0; i < end; i++)
            fact *= p*e - i - 1;

        // when picking the direction to reduce
        // we want to avoid reducing the factor coming from Frob(w)
        Vec<int64_t> shift;
        if( e > m)
            shift = w;
        else
            shift.SetLength(n + 1, int64_t(0));

        for(Hit = H.cbegin();  Hit != H.end(); Hit++)
        {
            const Vec<int64_t> &u = Hit->first;
            const Vec<R> &G = Hit->second;
            assert_print(u[0], ==, e);
            if( verbose > 2)
                cout <<"\t\treducing u = "<<u<<endl;

            for(i = 0; i < tuple_list[1].length(); i++)
                if( min_P(u - shift - tuple_list[1][i]) <= e - shift[0] - 1)
                    break;
            if( i >= tuple_list[1].length() )
            {
                cout << "couldn't find a v for u = "<<u<<endl;
                abort();
            }
            Vec<int64_t> &v = tuple_list[1][i];
            Vec<int64_t> dest = p*u - end*v;
            Vec<R> Gnew;
            R den;
            reduce_vector(Gnew, den, dest, v, end, G, method);
            assert_print(den, ==, 1);
            Vec<int64_t> final_dest;
            if( e > 1)
            {
                final_dest = u - v;
                assert_print(p*final_dest, ==, dest);
            }
            else
            {
                final_dest = u;
                assert_print(final_dest, ==, dest)
            }
            if( not IsZero(Gnew) )
            {
                Hnewit = Hnew.find(final_dest);
                if( Hnewit == Hnew.end())
                    Hnew[final_dest] = Gnew;
                else
                    Hnewit->second += Gnew;
            }
        }
        //end of loop in H
        H.swap(Hnew);

    }
    //end of loop in e


    if(F.length() != dim_dr_Y)
        F.SetLength(dim_dr_Y);
    for(i = 0; i < dim_dr_Y; i++)
            F[i] = R(0);

    assert_print(last_reduction_den, ==, 1);

    for(Hit = H.cbegin(); Hit != H.cend(); Hit++)
        F += last_reduction * (inclusion_matrices[ tuple_dict[1][Hit->first] ] * Hit->second);

    // F *= p^{n - 1}/ factorial(p * (m + N - 1) - 1)
    R fact_padic;
    int64_t val;
    factorial_padic(fact_padic, val, p * (m + N - 1) - 1, p);
    F *= inv(fact_padic);
    if(n - 1 >= val)
        F *= conv<R>(power_ZZ(p, n - 1 - val));
    else
    {
        R ppower = conv<R>(power_ZZ(p, val - (n - 1)));
        R q, r;
        for(i = 0; i < F.length(); i++)
        {

            divrem_lift(q, r, F[i], ppower);
            assert( IsZero(r) );
            F[i] = q;
        }
    }
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;

    if(verbose > 0)
    {
        cout <<"dr<R>::frob_monomial("<<w<<", "<<N<<") done!"<<endl;
        printf("Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );
    }

}


#endif
