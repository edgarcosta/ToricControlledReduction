// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

//see hypersurface.h for more details
void zeta_and_frob_Fp(Vec<ZZ> &zeta, Mat<ZZ> &F, Vec<int64_t> &hodge_numbers, Vec<int64_t>& r_vector, const int64_t &p, const map< Vec<int64_t>, ZZ, vi64less> &f, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const int64_t &verbose)
{
    int64_t a, weight, precision;
    Vec<int64_t> charpoly_prec, N;
    a = 1;
    
    
    //get the hodge numbers, to determine working precision
    {
        zz_pPush push(p);
        dr<zz_p> drminimal(p, conv<map< Vec<int64_t>, zz_p, vi64less> >(f), AP, bP, verbose, true);
        weight = drminimal.n - 1;
        hodge_numbers = drminimal.hodge_numbers;
    }

    //compute r_vector and charpoly_prec
    relative_precision(r_vector, charpoly_prec, hodge_numbers, weight, p, a); 
    //deduce N_vector
    N_vector(N, r_vector, weight, p);
    precision = working_precision(N, r_vector, p);
    if(verbose > 0)
    {
        print(hodge_numbers);
        print(r_vector);
        print(N);
        print(charpoly_prec);
        print(precision);
    }

    ZZ mod = power_ZZ(p, precision);

    if(NumBits(mod) < NTL_SP_NBITS)
    {
        zz_pPush push(conv<long>(mod));
        dr<zz_p> dr(p, conv<map< Vec<int64_t>, zz_p, vi64less> >(f), AP, bP, verbose, false);
        Mat<zz_p> Fp;
        dr.frob_matrix(Fp, N);
        F = conv< Mat<ZZ> >(Fp);
    }
    else
    {
        ZZ_pPush push(mod);
        dr<ZZ_p> dr(p, conv<map< Vec<int64_t>, ZZ_p, vi64less> >(f), AP, bP, verbose, false);
        Mat<ZZ_p> Fp;
        dr.frob_matrix(Fp, N);
        F = conv< Mat<ZZ> >(Fp);
    }
    if(verbose > 0)
        print(F);
    charpoly_frob(zeta, F, charpoly_prec, weight, p, a);
}


/*
 * same as above, but the input is given through string in the following format:
 *      p
 *      f.keys()
 *      f.values()
 *      AP
 *      bP
 */
void zeta_and_frob_Fp(Vec<ZZ> &zeta, Mat<ZZ> &F, Vec<int64_t> &hodge_numbers, Vec<int64_t>& r_vector, const char* input, const int64_t &verbose)
{
    int64_t p;
    map< Vec<int64_t>, ZZ, vi64less> f;
    Mat<int64_t> AP;
    Vec<int64_t> bP;
    stringstream buffer;
    buffer << input;
    buffer >> p;
    buffer >> f;
    buffer >> AP;
    buffer >> bP;
    zeta_and_frob_Fp(zeta, F, hodge_numbers, r_vector, p, f, AP, bP, verbose);
}


