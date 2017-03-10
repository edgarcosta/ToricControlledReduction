// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

//see hypersurface.h for more details


void frob_Fq(Mat<ZZX> &F, Vec<int64_t> &hodge_numbers, Vec<int64_t>& r_vector, Vec<int64_t> &charpoly_prec, const int64_t p, const ZZX &fE, map< Vec<int64_t>, ZZX, vi64less> &f, map< Vec<int64_t>, ZZX, vi64less> &ffrob, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const bool &verbose)
{
    int64_t a, weight, precision;
    Vec<int64_t> N;
    a = deg(fE);
    
    
    //get the hodge numbers, to determine working precision
    {
        zz_pPush push1(p);
        {
            zz_pEPush push2(conv<zz_pX>(fE));
            dr<zz_pE> drminimal(p, conv<map< Vec<int64_t>, zz_pE, vi64less> >(f), conv<map< Vec<int64_t>, zz_pE, vi64less> >(ffrob), AP, bP, verbose, true);
            weight = drminimal.n - 1;
            hodge_numbers = drminimal.hodge_numbers;
        }
    }

    //compute r_vector and charpoly_prec
    relative_precision(r_vector, charpoly_prec, hodge_numbers, weight, p, a); 
    //deduce N_vector
    N_vector(N, r_vector, weight, p);
    if(verbose > 0)
    {
        print(hodge_numbers);
        print(r_vector);
        print(N);
        print(charpoly_prec);
    }
    precision = working_precision(N, r_vector, p);
    ZZ mod = power_ZZ(p, precision);
    if(NumBits(mod) < NTL_SP_NBITS)
    {
        zz_pPush push1(conv<long>(mod));
        {
            zz_pEPush push2(conv<zz_pX>(fE));
            dr<zz_pE> dr(p, conv<map< Vec<int64_t>, zz_pE, vi64less> >(f), conv<map< Vec<int64_t>, zz_pE, vi64less> >(ffrob), AP, bP, verbose, false);
            Mat<zz_pE> Fp;
            dr.frob_matrix(Fp, N);
            F = conv< Mat<ZZX> >(Fp);
        }
    }
    else
    {
        ZZ_pPush push1(mod);
        {
            ZZ_pEPush push2(conv<ZZ_pX>(fE));
            dr<ZZ_pE> dr(p, conv<map< Vec<int64_t>, ZZ_pE, vi64less> >(f), conv<map< Vec<int64_t>, ZZ_pE, vi64less> >(ffrob), AP, bP, verbose, false);
            Mat<ZZ_pE> Fp;
            dr.frob_matrix(Fp, N);
            F = conv< Mat<ZZX> >(Fp);
        }
    }
    if(verbose > 0)
        print(F);
}
