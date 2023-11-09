// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "linear_algebra.h"

ZZ trace(const Mat<ZZ> &M)
{
    int64_t i;
    ZZ result;
    assert_print( M.NumRows(), ==, M.NumCols() );
    for(i = 0; i < M.NumRows(); i++)
        result += M[i][i];

    return result;
}
