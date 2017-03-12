// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

int main()
{
    bool val = true;
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;

// elliptic curve 11.a p = 17
const char ecp[] = "17 \n[ [0 2][0 0][3 0][2 0][0 1] ]\n[1 12 16 1 1] \n[ [ 0  1][ 1  0][-2 -3] ] \n[0 0 6] \n";
const char ecp_zeta[] = "[17 2 1]";
cout << "Testing: elliptic curve 11.a p = 17" <<endl;
user_time = get_cpu_time();
get_timestamp(&wtime1);
if( test_Fp(ecp, ecp_zeta, 0))
cout <<"PASS"<<endl;
else
val = false;
get_timestamp(&wtime2);
wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
user_time = get_cpu_time() - user_time;
printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

// genus 2 curve 11664.a  p = 31
const char g2p[] = "31 \n[ [6 0][0 0][0 2] ]\n[1 23 1] \n[ [ 1  0][ 0  1][-1 -3] ] \n[0 0 6] \n";
const char g2p_zeta[] = "[961 0 46 0 1]";
cout << "Testing: genus 2 curve 11664.a  p = 31" <<endl;
user_time = get_cpu_time();
get_timestamp(&wtime1);
if( test_Fp(g2p, g2p_zeta, 0))
cout <<"PASS"<<endl;
else
val = false;
get_timestamp(&wtime2);
wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
user_time = get_cpu_time() - user_time;
printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

// C_{3,4} curve p = 43
const char c34[] = "43 \n[ [0 3][0 0][4 0][1 0] ]\n[42 1 1 1] \n[ [ 0  1][ 1  0][-3 -4] ] \n[0 0 12] \n";
const char c34_zeta[] = "[79507 27735 6579 1258 153 15 1]";
cout << "Testing: C_{3,4} curve p = 43" <<endl;
user_time = get_cpu_time();
get_timestamp(&wtime1);
if( test_Fp(c34, c34_zeta, 0))
cout <<"PASS"<<endl;
else
val = false;
get_timestamp(&wtime2);
wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
user_time = get_cpu_time() - user_time;
printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

// plane curve, genus 3 p = 17
const char pg3[] = "17 \n[ [0 3][0 0][3 1][2 2][1 3][1 2][0 1][4 0][1 1][0 2][2 0][0 4] ]\n[7 1 14 16 1 16 6 15 15 15 5 1] \n[ [-1 -1][ 0  1][ 1  0] ] \n[4 0 0] \n";
const char pg3_zeta[] = "[4913 867 221 27 13 3 1]";
cout << "Testing: plane curve, genus 3 p = 17" <<endl;
user_time = get_cpu_time();
get_timestamp(&wtime1);
if( test_Fp(pg3, pg3_zeta, 0))
cout <<"PASS"<<endl;
else
val = false;
get_timestamp(&wtime2);
wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
user_time = get_cpu_time() - user_time;
printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

// genus 9 curve p = 37
const char g9[] = "37 \n[ [-1 -1][ 0  4][ 3  0] ]\n[1 1 36] \n[ [ 5 -1][-1  4][-4 -3] ] \n[4 3 12] \n";
const char g9_zeta[] = "[129961739795077 0 31612315085289 0 3417547576788 0 215521018356 0 8737338582 0 236144286 0 4254852 0 49284 0 333 0 1]";
cout << "Testing: genus 9 curve p = 37" <<endl;
user_time = get_cpu_time();
get_timestamp(&wtime1);
if( test_Fp(g9, g9_zeta, 0))
cout <<"PASS"<<endl;
else
val = false;
get_timestamp(&wtime2);
wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
user_time = get_cpu_time() - user_time;
printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

// hyperelliptic curve of genus 10 p = 13
const char g10[] = "13 \n[ [22  0][ 1  0][ 0  0][ 0  2] ]\n[1 1 1 12] \n[ [  0   1][  1   0][ -1 -11] ] \n[0 0 22] \n";
const char g10_zeta[] = "[137858491849 31813498119 8973037931 3074677333 323396203 24876631 10995985 12881011 3433911 999635 664001 76895 20319 5863 385 67 67 49 11 3 1]";
cout << "Testing: hyperelliptic curve of genus 10 p = 13" <<endl;
user_time = get_cpu_time();
get_timestamp(&wtime1);
if( test_Fp(g10, g10_zeta, 0))
cout <<"PASS"<<endl;
else
val = false;
get_timestamp(&wtime2);
wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
user_time = get_cpu_time() - user_time;
printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

// hyperelliptic curve of genus 10 p = 103
const char g10_103[] = "103 \n[ [22  0][ 1  0][ 0  0][ 0  2] ]\n[1 1 1 102] \n[ [  0   1][  1   0][ -1 -11] ] \n[0 0 22] \n";
const char g10_103_zeta[] = "[134391637934412192049 5219092735316978332 -152012409766513932 -105769152426538820 186272158258524 -534123936992982 124637044805304 8798421444054 969104502946 -108512162008 -4925241640 -1053516136 91347394 8051802 1107384 -46074 156 -860 -12 4 1]";
cout << "Testing: hyperelliptic curve of genus 10 p = 103" <<endl;
user_time = get_cpu_time();
get_timestamp(&wtime1);
if( test_Fp(g10_103, g10_103_zeta, 0))
cout <<"PASS"<<endl;
else
val = false;
get_timestamp(&wtime2);
wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
user_time = get_cpu_time() - user_time;
printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

return not val;
}
