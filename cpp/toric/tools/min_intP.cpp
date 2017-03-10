// Copyright 2017 Edgar Costa
// See LICENSE file for license details.


#include "tools.h"

//returns the minimal k such that v \in k*int(P)
int64_t min_intP(const Mat<int64_t> &AP, const Vec<int64_t> &bP, const Vec<int64_t> &v)
{
    //v \in k*int(P) iff AP*v + k*bP > 0;
    int64_t i, j, max, tmp;
    Vec<int64_t> w;
    w.SetLength(bP.length(), 0);
    // v[0] is just a dumb variable to potentially keep track the degree where v came from
    assert_print(AP.NumCols() + 1, ==, v.length());
    assert_print(AP.NumRows(), ==,  bP.length());
    // w = Ap * v[1:]
    for(i = 0; i < AP.NumRows(); i++)
        for(j = 0; j < AP.NumCols(); j++)
            w[i] += AP[i][j] * v[j + 1];
    
    max = -1;
    for(i = 0; i < bP.length(); i++)
    {
        if( (w[i] < 1) and (bP[i] == 0) )
            return INT64_MAX;//a practical +Infinity

        if( bP[i] != 0)
        {
            tmp = ceil(static_cast<double>( -(w[i] - 1) ) / bP[i] );
            if( tmp > max )
                max = tmp;
        }
    }
    return max;
}
