// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"
#include <cmath> // ceil, log
/*
 * Input:
 *  - hodge_numbers, hodge_numbers[i] = dim PH^(weight - i, i)(X) 
 *  - weight, weight of motive = dim(X) = n - 1
 *  - p, the characteristic of the field
 *  - a, q = p^a
 *
 * Output:
 *  - r_vector, the number of digits above the Hodge polygon necessary to deduce the correct characteristic polynomial, where r_vector[i] corresponds to PH^(weight - i, i)(X)
 *  - charpoly_prec, if compute the characteristic polynomial of a matrix representing an approximation of the Frobenius acting on PH^weight (X) with enough digits (according to r_vector), then the ith coefficient will be correct mod p^charpoly_prec[i]
 */
void relative_precision(Vec<int64_t> &r_vector, Vec<int64_t> &charpoly_prec, const Vec<int64_t> &hodge_numbers, const int64_t &weight, const int64_t &p, const int64_t &a)
{
    r_vector.kill();
    r_vector.SetLength( hodge_numbers.length(), int64_t(0));
    Vec<int64_t> HP, slope;
    hodge_polygon(HP, slope, hodge_numbers);
    charpoly_prec.kill();
    charpoly_prec.SetLength(HP.length(),int64_t(0));
    int64_t Pdeg, max_digits, i, j;
    Pdeg = HP.length() - 1;
    max_digits = 0;

    for(i = 1; 2*i < Pdeg + 2; i++)
    {
        int64_t r = ceil( log( double(2 * Pdeg)/i )/log((double)p) + a * i * weight * 0.5) - HP[i];
        if( r >= max_digits )
        {
            max_digits = r;
            for(j = 0; j < slope[i] + 1; j++)
                r_vector[j] = r;
            for(j = slope[i] + 1; j < hodge_numbers.length(); j++)
            {
                r--;
                r_vector[j] = r;
            }
        }
    }
    reverse(r_vector.begin(), r_vector.end());
    for(i = 0; i < HP.length(); i++)
        charpoly_prec[i] = HP[i] + max_digits;
    reverse(charpoly_prec.begin(), charpoly_prec.end());
}
