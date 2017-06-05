// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

int main()
{
    bool val = true;
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
    const char cubic[] = "47 \n[ [0 3 0][0 1 2][1 1 0][0 0 0][0 0 1][2 0 0][1 0 2][1 1 1][1 2 0][3 0 0][1 0 0][0 2 1][2 0 1][0 0 3][1 0 1][0 1 1][0 0 2][2 1 0][0 1 0][0 2 0] ]\n[42 26 5 45 26 26 26 5 26 41 26 26 26 43 5 5 26 26 26 26] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[3 0 0 0] \n";
    const char cubic_zeta[] = "[ -10779215329 0 0 0 0 0 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: cubic surface p = 47" <<endl;
    if( test_Fp(cubic, cubic_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );


    const char cubic2[] = "1009 \n[ [0 3 0][0 1 2][1 1 0][0 0 0][0 0 1][2 0 0][1 0 2][1 1 1][1 2 0][3 0 0][1 0 0][0 2 1][2 0 1][0 0 3][1 0 1][0 1 1][0 0 2][2 1 0][0 1 0][0 2 0] ]\n[1004 988 967 1007 988 988 988 967 988 1003 988 988 988 1005 967 967 988 988 988 988] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[3 0 0 0] \n";
    const char cubic2_zeta[] = "[1055229678769825441 1045817322864049 0 -1027243729 0 1009 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: cubic surface p = 1009" <<endl;
    if( test_Fp(cubic2, cubic2_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );
    
    const char cubic3[] = "65537 \n[ [0 3 0][0 1 2][1 1 0][0 0 0][0 0 1][2 0 0][1 0 2][1 1 1][1 2 0][3 0 0][1 0 0][0 2 1][2 0 1][0 0 3][1 0 1][0 1 1][0 0 2][2 1 0][0 1 0][0 2 0] ]\n[65532 65516 65495 65535 65516 65516 65516 65495 65516 65531 65516 65516 65516 65533 65495 65495 65516 65516 65516 65516] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[3 0 0 0] \n";
    const char cubic3_zeta[] = "[79235416345888816038194577409 0 -18447869999386460161 0 -4295098369 0 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: cubic surface p = 2^16 + 1" <<endl;
    if( test_Fp(cubic3, cubic3_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

/*
    const char cubic4[] = "1048583 \n[ [0 3 0][0 1 2][1 1 0][0 0 0][0 0 1][2 0 0][1 0 2][1 1 1][1 2 0][3 0 0][1 0 0][0 2 1][2 0 1][0 0 3][1 0 1][0 1 1][0 0 2][2 1 0][0 1 0][0 2 0] ]\n[1048578 1048562 1048541 1048581 1048562 1048562 1048562 1048541 1048562 1048577 1048562 1048562 1048562 1048579 1048541 1048541 1048562 1048562 1048562 1048562] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[3 0 0 0] \n";
    const char cubic4_zeta[] = "[1329281237998693845036542985209236369 -1267692913196851222112644383143 0 0 0 -1048583 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: cubic surface p = 2^20 + 7" <<endl;
    if( test_Fp(cubic4, cubic4_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );
*/



    return not val;
}
