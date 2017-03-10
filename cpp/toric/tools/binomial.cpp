// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "tools.h"

ZZ binomial(int64_t n, int64_t k)
{
    assert(n>=k);
    ZZ num, den;
    int64_t i;
    num = 1; // n!/k! = (k+1) ... n
    den = 1; // (n-k)!
    for( i = k + 1; i <= n; i++)
    {
        num *= i;
    }
    for( i = 1; i <= n - k; i++)
    {
        den *=i;
    }
    return num/den;
}

