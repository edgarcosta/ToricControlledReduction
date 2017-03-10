// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "linear_algebra.h"

//computes the characteristic polynomial of M using Newton identities
// = det(I*X - M) = X^n - trace(M)*X^(n - 1) + ...  
Vec<ZZ> charpoly(const Mat<ZZ> &M)
{
    int64_t i, k, n;
    int64_t rem;
    Vec<ZZ> result, e, p;
    Mat<ZZ> powerM;
    ZZ sump, sumn;
    n = M.NumRows();
    assert( n == (int64_t) M.NumCols() );

    result.SetLength(n+1);
    e.SetLength(n+1);
    p.SetLength(n+1);
    set(e[0]);
    set(result[n]);
    p[1] = trace(M);

    powerM = M;

    for(i = 2; i <=n; i++)
    {
        powerM *= M; //powerM = M^i;
        p[i] = trace(powerM);
    }

    for( k = 1; k <= n; k++)
    {
        clear(sumn);
        clear(sump);
        for( i = 1; i <= k ; i+=2)
        {
            sump += e[k-i]*p[i];
        }
        for( i = 2; i <= k; i+=2)
        {
            sumn += e[k-i]*p[i];
        }

        rem = DivRem(e[k], sump - sumn, k);
        assert(rem == 0);
        result[n-k] = (k%2 == 0)? e[k]: -e[k];
    }
    return result;
}
