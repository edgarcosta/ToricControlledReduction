// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

bool test_Fq(const char *input, const char *output, const int64_t &verbose)
{
    stringstream buffer;
    
    Mat<ZZX> F, F_expected;
    Vec<int64_t> hodge_numbers, r_vector, charpoly_prec;
    frob_Fq(F, hodge_numbers, r_vector, charpoly_prec, input, verbose);

    buffer << output;
    buffer >> F_expected;
    
    if(F_expected != F)
    {
        cout << "FAIL" <<endl;
        print(F);
        print(F_expected);
        return false;
    }
    return true;
}

