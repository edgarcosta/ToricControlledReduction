// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include <set>
#include "tools.h"



void integral_points(Vec< Vec< Vec< int64_t > > > &tuple_list, Vec< Vec< Vec<int64_t> > > &tuple_int_list, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const Vec< Vec<int64_t> > &Pvertices, const int64_t &N)
{
    /*
     * Assumes that P is the convex hull of Pvertices and 0
     * and P is a normal Polytope, ie, we can deduce k*P = P_1 + ..... + P_1
     */
    
    /* 
     * First we compute a cube, prod( [Ai, Bi]), such that P < prod( [Ai, Bi])
     * From it we will deduce the integral points of P, and from there we use normality to deduce k*P
     */
    assert(N > 0);
    Vec<int64_t> A, B, zA, zB, v;
    int64_t i, j, n, k;
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



    zA.SetLength(1, 1);
    zA.append(A);

    // NB = [0] + N * B
    zB.SetLength(1, 1);
    zB.append(B);

    v = zA;

    std::set< Vec< int64_t >, vi64less> S;
    std::set< Vec< int64_t >, vi64less>::const_iterator Sit;
    // v loop over the box  prod(  [-Ai, Bi] )
    while( true )
    {
        for( ;v[n] <= zB[n]; v[n]++)
        {
            //do stuff
            if( min_P(AP, bP, v) <= 1 )
                S.insert(v);
            
        }
        v[n] = zB[n];

        if(v == zB)
            break;
        else
        {
            // v[0] is a dummy variable
            for(i = n; i >= 1; i--)
            {
                if( v[i] == zB[i] )
                    v[i] = zA[i];
                else
                {
                    v[i]++;
                    break;
                }
            }
        }
    }

    tuple_list.SetLength(N);
    tuple_int_list.SetLength(N);
    //allocate enough space
    tuple_list[0].SetLength(1);
    tuple_list[0][0].SetLength(n + 1, 0);
    tuple_list[1].SetLength(S.size());
    tuple_int_list[0].SetLength(0);
    tuple_int_list[1].SetMaxLength(S.size());

    for(i = 0, Sit = S.cbegin(); Sit != S.cend(); Sit++, i++)
        tuple_list[1][i] = *Sit;
    

    Vec< Vec<int64_t> > buffer;
    for(k = 2; k < N; k++)
    {
        S.clear();
        for(i = 0; i < tuple_list[1].length(); i++)
        {
            v = tuple_list[1][i];
            buffer.SetLength(tuple_list[k-1].length());
            for(j = 0; j < tuple_list[k-1].length(); j++)
                buffer[j] =  v + tuple_list[k-1][j];
            S.insert(buffer.begin(), buffer.end());
        }
        tuple_list[k].SetLength(S.size());
        tuple_list[k].FixAtCurrentLength();
        for(i = 0, Sit = S.cbegin(); Sit != S.cend(); Sit++, i++)
            tuple_list[k][i] = *Sit;
        tuple_int_list[k].SetMaxLength(tuple_list[k].length());
    }
    
    for(i = 0; i < tuple_list[N - 1].length(); i++)
    {
        v = tuple_list[N-1][i];
        for(k = min_intP(AP, bP, v); k < N; k++)
        {
            v[0] = k;
            tuple_int_list[k].append(v);
        }
    }
    for(k = 1; k < N; k++)
    {
        tuple_int_list[k].SetMaxLength(tuple_int_list[k].length());
        tuple_int_list.FixAtCurrentLength();
    }
}
