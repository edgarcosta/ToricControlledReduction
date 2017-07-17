// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

bool test_Fp(const char *input, const char *output, const int64_t &verbose)
{
    stringstream buffer;
    
    Vec<ZZ> zeta, zeta_expected;
    Mat<ZZ> F;
    Vec<int64_t> hodge_numbers, r_vector;
    zeta_and_frob_Fp(zeta, F, hodge_numbers, r_vector, input, verbose);

    buffer << output;
    buffer >> zeta_expected;
    
    if(zeta_expected != zeta)
    {
        cout << "FAIL" <<endl;
        print(zeta);
        print(zeta_expected);
        print(hodge_numbers);
        return false;
    }
    return true;
}

