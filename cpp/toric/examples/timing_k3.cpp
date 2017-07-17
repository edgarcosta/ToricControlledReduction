// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

int main()
{
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
    
    int64_t p;
    int64_t precision = 4;
    Vec<int64_t> N;
    Vec<int64_t> ps;
    ps.SetLength(2);
    ps[0] = 1009;
    ps[1] = 10007;
    ps[0] = 109;
    ps[1] = 1009;
    N.SetLength(3, 3);
    N[2] = 2;

    
    stringstream buffer;
    const char input[] = "[ [0 3 0][3 1 0][0 0 0][0 0 3][0 2 2][0 2 1][0 1 0][4 0 0][1 0 2][0 1 1][2 0 0][0 1 2][1 1 0][1 3 0][2 2 0][2 0 1][1 0 1][1 0 0][2 0 2][1 2 1][3 0 0][1 0 3][0 3 1][0 4 0][0 0 1][1 1 1][0 1 3][0 0 2][3 0 1][1 2 0][2 1 0][2 1 1][0 0 4] ] \n [-7 -10 -5 1 6 3 3 -9 -1 8 -1 -3 -4 2 -9 9 -3 4 8 -2 8 2 3 -7 4 -9 7 -6 6 -8 1 9 9]\n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n";
    buffer << input;
    map< Vec<int64_t>, int64_t, vi64less> fZZ;
    buffer >> fZZ;

    Mat<int64_t> AP;
    Vec<int64_t> bP;
    buffer >> AP;
    buffer >> bP;

    
    int64_t verbose = 0;

    for(int64_t i = 0; i < ps.length(); i++)
    {
        p = ps[i];
        ZZ mod = power_ZZ(p, precision);
        zz_pPush push(conv<long>(mod));
        dr<zz_p> dr(p, conv<map< Vec<int64_t>, zz_p, vi64less> >(fZZ), AP, bP, verbose, false);
        for(int64_t method = 0; method < 3; method++)
        {
            user_time = get_cpu_time();
            get_timestamp(&wtime1);
            Vec<zz_p> res;
            dr.frob_monomial(res, dr.basis_dr_X[1][0], N[0], 0);
            get_timestamp(&wtime2);
            wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
            user_time = get_cpu_time() - user_time;
            printf("p = %lld\tmethod=%lld\tTime: CPU %.2f s, Wall: %.2f s\n", (long long) p, (long long) method, user_time, wall_time );
        }
        cout<<endl;
    }
}
