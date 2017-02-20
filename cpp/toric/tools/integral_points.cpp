// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "tools.h"



void integral_points(Vec< Vec< Vec< int64_t > > > &local_tuple_list, Vec< Vec< Vec<int64_t> > > &local_tuple_int_list, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const Vec< Vec<int64_t> > &Pvertices, const int64_t &N)
{
    //assumes that P is the convex hull of Pvertices and 0
    //hence we will use to Pvertices to figure out Ai and Bi such that P < prod( [Ai, Bi])
    //and then loop over all monomials in prod( (n+1) [Ai, Bi]) to initialize tuple_list[i] i <= n + 1;
    Vec<int64_t> A, B, NA, NB, v;
    int64_t i, j;
    int64_t n, Pmax_volume, vol, k, kint;
    n = Pvertices[0].length();
    A.SetLength(n, 0);
    B.SetLength(n, 0);
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
    for(i = 0; i < n; i++)
        assert_print(A[i], < , B[i]);

    //w = B - A + (1,...,1);
    sub(v, B, A);
    for(i = 0; i < n; i++)
        v[i]++;
    // prod(w) = Vol( prod( [-Ai, Bi]) )
    Pmax_volume = prod(v);

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
    
    // NA = [0] + N * A
    mul(v, N, A);
    NA.SetLength(1, 0);
    NA.append(v);
        
    // NB = [0] + N * B
    mul(v, N, B);
    NB.SetLength(1, 0);
    NB.append(v);

    
    v = NA;

    // v loop over the box  prod( (n+1) [-Ai, Bi] )
    while( true )
    {
        for( ;v[n] <= NB[n]; v[n]++)
        {
            //do stuff
            
            k = min_P(AP, bP, v);
            kint = min_intP(AP, bP, v);
            for(j = k; j < N; j++)
            {
                v[0] = j;
                local_tuple_list[j].append(v);
                if( j >= kint )
                    local_tuple_int_list[j].append(v);
            }
        }
        v[n] = NB[n];
        v[0] = 0;

        if(v == NB)
            break;
        else
        {
            // v[0] is a dummy variable
            for(i = n; i >= 1; i--)
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
