// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "tools.h"



void integral_points(Vec< Vec< int64_t > > local_tuple_list, Vec< Vec<int64_t> > local_tuple_int_list, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const Vec<Vec<int64_t>> &Pvertices, const int64_t N)
{
    //assumes that P is the convex hull of Pvertices and 0
    //hence we will use to Pvertices to figure out Ai and Bi such that P < prod( [-Ai, Bi])
    //and then loop over all monomials in prod( (n+1) [-Ai, Bi]) to initialize tuple_list[i] i <= n + 1;
    Vec<int64_t> A, B, NA, NB, v, w;
    int64_t i, j;
    int64_t n, Pmax_volume, vol, k, kint;
    n = Pvertices[0].length();
    A.SetLength(n);
    B.SetLength(n);
    for(j = 0; j < Pvertices.length(); j++)
    {
        for(i = 0; i < n; i++)
        {
            if(Pvertices[j][i] < A[i])
                A[i] = Pvertices[j][i];

            if(Pvertices[j][i] > B[i])
                B[i] = Pvertices[j][i];
        }
    }
    //w = B - A + (1,...,1);
    sub(w, B, A);
    for(i = 0; i < n; i++)
        w[i]++;
    // prod(w) = Vol( prod( [-Ai, Bi]) )
    Pmax_volume = prod(w);

    local_tuple_list.SetLength(N);
    local_tuple_int_list.SetLength(N);
    //allocate enough space
    local_tuple_list[0].SetMaxLength(1);
    local_tuple_int_list[0].SetMaxLength(0);
    for(i =1; i < N; i++)
    {
        vol = Pmax_volume * pow(i, n);
        local_tuple_list[i].SetMaxLength(vol);
        local_tuple_int_list[i].SetMaxLength(vol);
    }
    
    // Ab = (n + 1) * A
    mul(NA, n + 1, A);
    // Bn1 = (n + 1) * B
    mul(NB, n + 1, B);

    v = NA;
    w.SetLength(n + 1);
    // v loop over the box  prod( (n+1) [-Ai, Bi] )
    while( true )
    {
        for( ;v[n-1] <= NB[n-1]; v[n-1]++)
        {
            //do stuff
            for(j = 1; j < n + 1; j++)
                w[j] = v[j];
            
            k = min_P(AP, bP, w);
            kint = min_intP(AP, bP, w);

            //k or kint = -1, means that w \notin k*P for any k
            if( k == -1 )
                k = n + 2;
            if( kint == -1 )
                k = n + 2;

            for(j = k; j < n + 2; j++)
            {
                w[0] = j;
                local_tuple_list.append(w);
                if( j > kint )
                    local_tuple_list.append(w);
            }

        }
        v[n-1] = B[n-1];

        if(v == B)
            break;
        else
        {
            for(i = n - 1; i >= 0; i--)
            {
                if( v[i] == NB[i] )
                    v[i] = NA[i];
                else
                {
                    v[i]++;
                    break;
                }
            }
        }
    }
}
