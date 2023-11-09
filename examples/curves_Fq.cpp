// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

int main()
{
    bool val = true;
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;


    // elliptic curve 11.a
    const char ecq[] = "17 \n[3 16 1] \n[ [0 2][0 0][3 0][2 0][0 1] ]\n[ [1] [12] [16] [1] [1]  ] \n[ [2 0][0 0][3 0][0 2][0 1] ]\n[ [1] [12] [16] [1] [1]  ] \n[ [ 0  1][ 1  0][-2 -3] ] \n[0 0 6] \n";
    const char ecq_F[] = "[ [ [4760] [195] ] [ [2023] [151] ] ] ";
    cout << "Testing: elliptic curve 11.a" <<endl;
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    if( test_Fq(ecq, ecq_F, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );


    // genus 2 curve 11664.a
    const char g2q[] = "31 \n[3 29 1] \n[ [6 0][0 0][0 2] ]\n[ [1] [23] [1]  ] \n[ [6 0][0 0][0 2] ]\n[ [1] [23] [1]  ] \n[ [ 1  0][ 0  1][-1 -3] ] \n[0 0 6] \n";
    const char g2q_F[] = "[ [ [276055] [] [] [] ] [ [] [647466] [] [] ] [ [] [] [7932] [] ] [ [] [] [] [21859] ] ] ";
    cout << "Testing: genus 2 curve 11664.a" <<endl;
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    if( test_Fq(g2q, g2q_F, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );
    /*
     * Same curves under a change of variables:
     * xq = xq * (1 - Fq.gen())
     * yq = (1 + Fq.gen())*yq
     */
    // elliptic curve 11.a
    const char ecq1[] = "17 \n[3 16 1] \n[ [0 2][0 0][3 0][2 0][0 1] ]\n[ [1] [12] [16] [1] [1]  ] \n[ [2 0][0 0][3 0][0 2][0 1] ]\n[ [1] [12] [16] [1] [1]  ] \n[ [ 0  1][ 1  0][-2 -3] ] \n[0 0 6] \n";
    const char ecq1_F[] = "[ [ [4760] [195] ] [ [2023] [151] ] ] ";
    cout << "Testing: elliptic curve 11.a" <<endl;
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    if( test_Fq(ecq1, ecq1_F, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );


    // genus 2 curve 11664.a
    const char g2q1[] = "31 \n[3 29 1] \n[ [6 0][0 0][0 2] ]\n[ [1] [23] [1]  ] \n[ [6 0][0 0][0 2] ]\n[ [1] [23] [1]  ] \n[ [ 1  0][ 0  1][-1 -3] ] \n[0 0 6] \n";
    const char g2q1_F[] = "[ [ [276055] [] [] [] ] [ [] [647466] [] [] ] [ [] [] [7932] [] ] [ [] [] [] [21859] ] ] ";
    cout << "Testing: genus 2 curve 11664.a" <<endl;
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    if( test_Fq(g2q1, g2q1_F, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );


    return not val;
}
