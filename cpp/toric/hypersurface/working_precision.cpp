// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "tools.h"


// outputs =  max([ self.r_vector[m - 1] + max(factorial_p_adic(self.p*(m + self.N[m - 1] - 1) - 1, self.p)[0] - m + 1, 0) for m in range(1, self.n + 1) ])
int64_t working_precision(const Vec<int64_t> &N, const Vec<int64_t> &r_vector, const int64_t &p)
{
    Vec<int64_t> v;
    v.SetLength(r_vector.length());
    for(int64_t k = 0; k <  r_vector.length(); k++)
    {
        //k = m - 1
        //self.r_vector[m - 1] + max(factorial_p_adic(self.p*(m - 1 + self.N[m - 1]) - 1, self.p)[0] - (m - 1), 0)
        v[k] = r_vector[k] + max(int64_t(0), valuation_of_factorial(p * (k + N[k]) - 1, p) - k);
    }
    return max(v);
}
