// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"
#include <cmath>

void charpoly_frob(Vec<ZZ> &cp, const Mat<ZZ> &frob_matrix, const Vec<int64_t> &charpoly_prec, const int64_t &weight, const int64_t &p, const int64_t &a)
{
    int64_t i, k;
    cp = charpoly(frob_matrix);

    assert_print(cp.length(), ==, charpoly_prec.length());
    assert_print(cp[cp.length() - 1], ==, 1);

    int64_t sign, degree;
    degree = cp.length() - 1;
    int64_t halfdegree = ceil(degree*0.5) + 1;
    sign = 0;
    Vec<ZZ> mod;
    mod.SetLength( charpoly_prec.length() );
    for(i = 0; i <= degree; i++)
    {
        mod[i] =power_ZZ(p, charpoly_prec[i]);
        cp[i] = cp[i] % mod[i];
    }
    
    // figure out the sign
    if( weight % 2 == 1)
        // for odd weight the sign is always 1
        // it's the charpoly of a USp matrix
        // and charpoly of a symplectic matrix is reciprocal
        sign = 1;
    else
    {
        for(i = 0; i < degree + 1; i++)
        {
            ZZ p_power = power_ZZ(p, min(charpoly_prec[i], charpoly_prec[degree - i] + a * (degree - 2 * i)) * (weight /2) );
            if( cp[i] % p_power !=0 && cp[degree-i] % p_power != 0)
            {
                if (0 == (cp[i] + cp[degree - i] * power_ZZ(p, a * (degree - 2 * i)) * (weight /2) ) %  p_power )
                    sign = -1;
                else
                    sign = 1;

                assert_print(0, ==, (sign * cp[i] - cp[degree - i] * power_ZZ(p, a * (degree - 2 * i)) * (weight /2) ) %  p_power );
                break;
            }
        }
    }
    assert_print(sign, != , 0);
    cp[0] = sign * power_ZZ(p,  (a * degree * weight)/2);
    //calculate the i-th power sum of the roots and correct cp allong the way
    Vec<ZZ> e;
    e.SetLength(halfdegree);
    // e[k] = \sum x_{i_1} x_{i_2} ... x_{i_k} # where x_* are eigenvalues
    // and i_1 < i_2 ... < i_k
    for(k = 0; k < halfdegree; k++)
    {
        e[k] = (k%2 == 0)? cp[degree - k] : -cp[degree - k];
        if(k > 0)
            //verify if p^charpoly_prec[degree - k] > 2*degree/k * q^(w*k/2
            assert_print( log(k)/log(p) + charpoly_prec[degree-k],  > ,log(2*degree)/log(p) + a*0.5*weight*k );
    }
    //s[k] = \sum x_i ^k for i>0
    ZZ sum, pN;
    Vec<ZZ> s;
    s.SetLength(halfdegree);
    for(k = 1; k < halfdegree; k++)
    {
        /*
        * assume that s[i] and e[i] are correct for i < k
        * e[k] correct modulo mod[degree - k]
        * S = - sum (-1)^i e[k-i] * s[i] = sum (-1)^(i+1) e[k-i] * s[i]
        * s[k] = (-1)^k (S - k*e[k] ) 
        * ==> k*e[k] = (-1)^(k+1) s[k] + S 
        */
        sum = 0;
        for(i = 1; i < k ; i++)
            sum += ((i%2 == 1)? 1 : -1) * e[k-i] * s[i];
        s[k] = ( (k%2 == 0)? 1 : -1 ) * (sum - k * e[k]);
        //hence s[k] is correct modulo k*mod[degree - k] 
        pN = k * mod[degree - k];
        s[k] = s[k] % pN;
        // |x_i| = p^(w*0.5)
        // => s[k] <= degree*p^(a*w*k*0.5)
        // recall, 2*degree*p^(a*w*k*0.5) /k < mod[degree - k] 
        if( sqr(s[k]) > degree * degree * power_ZZ(p, a * weight * k) )
            s[k] = - ( (-s[k]) % pN );

        //now correct e[k] with:
        // (-1)^(k+1) s[k] + S = k*e[k] 
        sum = sum + ( (k%2 == 0)? -s[k] : s[k] );
        int64_t rem = DivRem(e[k], sum , k);
        assert_print(rem, ==, 0);
        cp[degree - k] = (k%2 == 0)? e[k]: -e[k];
        //either weight is even, or weight is odd and degree must be even
        if( degree >= 2*k )
            cp[k] =  sign * cp[degree - k] * power_ZZ(p,  (a* (degree - 2 * k) * weight)/2 );
    }
}

