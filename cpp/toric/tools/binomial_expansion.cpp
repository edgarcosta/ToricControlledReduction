// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "tools.h"

void binomial_expansion(Vec<ZZ> &v, const int64_t &a, const int64_t &b, const int64_t &n)
{
    v.kill();
    if(a == 0 or b == 0 or n == 0)
    {
        if( n == 0)
        {
            v.SetLength(1);
            v[0] = 1;
            return;
        }
        if( a == 0 and b == 0)
        {
            v.SetLength(1);
            v[0] = 0;
            return;
        }

        v.SetLength(n + 1, ZZ(0));
        if( b == 0)
        {
            power(v[0], a, n);
            return;
        }
        if( a == 0)
        {
            power(v[n], b, n);
            return;
        }
    }
    //compute the binomial coefficients
    v.SetLength(n + 1, ZZ(1));
    ZZ c(1);
    for(int64_t k = 1; 2*k < n + 1; k++)
    {
        c = (c * ( n - k + 1))/k;
        v[k] = c;
        v[n - k] = c;
    }
    
    ZZ apower(1);
    ZZ bpower(1);
    for(int64_t k = 1; k < n + 1; k++)
    {
        //apower = a^k
        //bpower = b^k
        apower *= a;
        bpower *= b;

        //v[k] = (n choose k) * a^(n-k) * b^k
        //v[n - k]  = (n choose k) * a^k * b^(n -k)
        v[k] *= bpower;
        v[n - k] *= apower;
    }
}

