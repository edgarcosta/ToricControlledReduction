// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "solve_system.h"
#include <NTL/LLL.h>

template<> void pivot_columns(Vec<int64_t> &res, const Mat<ZZ> &T)
{
    while(true)
    {
        //pick a random prime p ~ 2^23, so it can possibly use AVX and FMA instructions
        zz_pPush push(GenPrime_long(23));
        pivot_columns(res, conv< Mat<zz_p> >(T) );
        int64_t rank = res.length();
        Mat<ZZ> S;
        // S = T on top of rows forcing the nonpivot columns of T to be pivots of S
        S.SetDims(T.NumRows() + T.NumCols() - rank, T.NumCols());
        int64_t i, j, position;
        for(i = 0; i < T.NumRows(); i++)
            for(j = 0; j < T.NumCols(); j++)
                S[i][j] = T[i][j];
        //force the non pivot columns to be pivot
        //by adding the respective rows to S
        position = T.NumRows();
        i = 0;
        for(j = 0; j < rank; j++)
        {
            while(i < res[j] and i < T.NumCols())
            {
                S[position][i] = 1;
                i++;
                position++;
            }
            i++;
        }
        while(i < T.NumCols())
        {
                S[position][i] = 1;
                i++;
                position++;
        }
        assert( position == S.NumRows() );

        ZZ det2;
        int64_t rankS;
        rankS = image(det2, S);
        //if S is a full rank matrix, then the pivot matrices over Fp are indeed global pivots
        if(rankS == T.NumCols())
            break;
    }
}
