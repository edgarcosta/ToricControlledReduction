// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

// 11.a elliptic curve
// http://www.lmfdb.org/EllipticCurve/Q/11/a/1
int main()
{
    const char ecp[] = "17 \n[ [0 2][0 0][3 0][2 0][0 1] ]\n[1 12 16 1 1] \n[ [ 0  1][ 1  0][-2 -3] ] \n[0 0 6] \n";
    const char zexpected[] = "[17 2 1]";
    stringstream buffer;
    
    Vec<ZZ> zeta, zeta_expected;
    Mat<ZZ> F;
    Vec<int64_t> hodge_numbers, r_vector;
    zeta_and_frob_Fp(zeta, F, hodge_numbers, r_vector, ecp);

    buffer << zexpected;
    buffer >> zeta_expected;
    
    assert_print(zeta_expected, ==, zeta);
    return 0;
}
