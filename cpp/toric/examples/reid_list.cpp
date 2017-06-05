// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

int main()
{
    timestamp_type wtime1, wtime2;
    Vec<ZZ> zeta;
    Mat<ZZ> F;
    Vec<int64_t> hodge_numbers, r_vector;
    double wall_time, user_time;
    const char reid62_109[] = "109 \n[ [0 0 0][2 0 0][0 0 4][0 1 5][1 0 4][0 4 0] ]\n[1 1 1 1 1 1] \n[ [ 0  1  0][ 1  0  0][ 0  1 -1][-8 -5 -3][-4 -2 -1][ 0  0  1] ] \n[0 0 4 20 8 0] \n";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    zeta_and_frob_Fp(zeta, F, hodge_numbers, r_vector, reid62_109, 0);
    print(zeta);
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );
    
    const char reid62_10007[] = "10007 \n[ [0 0 0][2 0 0][0 0 4][0 1 5][1 0 4][0 4 0] ]\n[1 1 1 1 1 1] \n[ [ 0  1  0][ 1  0  0][ 0  1 -1][-8 -5 -3][-4 -2 -1][ 0  0  1] ] \n[0 0 4 20 8 0] \n";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    zeta_and_frob_Fp(zeta, F, hodge_numbers, r_vector, reid62_10007, 0);
    print(zeta);
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );


    const char reid62_49999[] = "49999 \n[ [0 0 0][2 0 0][0 0 4][0 1 5][1 0 4][0 4 0] ]\n[1 1 1 1 1 1] \n[ [ 0  1  0][ 1  0  0][ 0  1 -1][-8 -5 -3][-4 -2 -1][ 0  0  1] ] \n[0 0 4 20 8 0] \n";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    zeta_and_frob_Fp(zeta, F, hodge_numbers, r_vector, reid62_49999, 0);
    print(zeta);
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    return 0;
}
