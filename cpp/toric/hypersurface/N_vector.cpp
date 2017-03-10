// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

void N_vector(Vec<int64_t> &N, const Vec<int64_t> &r_vector, const int64_t &weight, const int64_t &p)
{
    if( p < 2*(weight + 1) + max(r_vector) )
    {
        cout <<"p is too small, only implemented for p >= 2*(weight + 1) + r"<<endl;
        print(p);
        print(2*(weight + 1));
        print(max(r_vector));
        cout<<"bye bye"<<endl;
        abort();
    }
    else
    {
        N.SetLength(r_vector.length());
        for(int64_t i = 0; i < N.length(); i++)
            N[i] = (weight + 1) + r_vector[i] - (i + 1);
    }
}
