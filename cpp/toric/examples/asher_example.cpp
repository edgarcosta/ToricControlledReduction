// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

int main()
{
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
    

    Vec<ZZ> zeta;
    Mat<ZZ> F;
    Vec<int64_t> hodge_numbers, r_vector, N, charpoly_prec;

    /*
     #asher's example
     P5 = Polyhedron([[0,0,0,0,0],[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
R4.<y_1, y_2, y_3, y_4, y_5> = QQ[]
R5.<y_0> = R4[]


AA = y_0^2*y_1+y_0*y_1^2+y_1^2*y_2+y_0*y_2^2+4*y_0^2*y_3+y_1^2*y_3+8*y_0*y_2*y_3+2*y_1*y_2*y_3+2*y_2^2*y_3+4*y_0*y_3^2+y_1*y_3^2+y_3^3+8*y_0*y_1*y_4+y_1^2*y_4+4*y_1*y_2*y_4+y_2^2*y_4+8*y_0*y_3*y_4+2*y_2*y_3*y_4+8*y_0*y_4^2+y_1*y_4^2+2*y_3*y_4^2+y_4^3+2*y_0^2*y_5+y_1^2*y_5+y_1*y_2*y_5+y_2^2*y_5+8*y_0*y_3*y_5+y_1*y_3*y_5+y_3^2*y_5+4*y_0*y_4*y_5+3*y_3*y_4*y_5+2*y_0*y_5^2+y_4*y_5^2
change_of_variables = (-y_1 - y_3 + y_4 + 3*y_5, -y_0 + 2*y_2 + y_3, -y_2 + y_3 + y_4 - 5*y_5, -2*y_1 + y_2 + y_3 + y_4 + 4*y_5, 9*y_0 + y_1 + y_3 - y_4 - y_5, -y_1 - y_2 + 5*y_3 + 4*y_4 - y_5)

AAnew = R4( AA(tuple(change_of_variables)).subs(y_0 = 1))
P = Polyhedron(AAnew.dict().keys())
if P == 3*P5:
    print repr(input_to_NTL_Fp(AAnew.change_ring(GF(47)).dict()))
*/
    const char input[] = "23 \n[ [2 0 0 0 0][0 0 0 0 0][0 3 0 0 0][0 0 2 0 0][0 1 1 1 0][0 2 1 0 0][0 0 1 2 0][0 0 2 1 0][0 1 1 0 0][0 0 0 2 1][1 0 2 0 0][1 0 0 1 1][0 0 0 2 0][1 1 0 1 0][0 0 0 1 2][0 0 0 3 0][0 0 1 0 1][1 1 0 0 0][1 1 0 0 1][1 0 1 0 1][0 1 2 0 0][0 1 0 1 0][1 0 0 1 0][0 0 1 0 0][0 2 0 0 0][0 1 0 0 1][1 0 0 2 0][0 0 0 0 1][1 0 0 0 0][1 0 1 1 0][3 0 0 0 0][0 0 3 0 0][0 1 0 0 2][0 1 0 2 0][2 0 0 1 0][2 0 1 0 0][1 0 0 0 1][0 0 0 1 0][0 0 1 1 1][0 0 1 0 2][1 0 1 0 0][0 1 0 1 1][1 0 0 0 2][0 1 1 0 1][0 0 0 1 1][0 1 0 0 0][2 1 0 0 0][1 1 1 0 0][1 2 0 0 0][0 2 0 0 1][0 2 0 1 0][0 0 0 0 3][0 0 1 1 0] ]\n[14 17 2 13 19 12 6 18 17 21 16 10 11 19 13 22 1 11 19 13 7 3 7 2 2 8 6 12 16 18 1 1 17 15 1 3 19 11 7 22 6 7 3 21 10 10 14 3 10 19 9 7 19] \n[ [ 0  0  0  0  1][-1 -1 -1 -1 -1][ 0  0  0  1  0][ 0  0  1  0  0][ 0  1  0  0  0][ 1  0  0  0  0] ] \n[0 3 0 0 0 0] \n";
    zeta_and_frob_Fp(zeta, F, hodge_numbers, r_vector, input, 100);
    hodge_numbers.SetLength(5, int64_t(0));
    hodge_numbers[0] = hodge_numbers[4] = 0;
    hodge_numbers[1] = hodge_numbers[3] = 1;
    hodge_numbers[2] = 20;

    N.SetLength(5, 0);
    N[0] = N[1] = N[2] = N[3] = N[4] = 4;
    r_vector.SetLength(5, int64_t(0));
    r_vector[0] = 0;
    r_vector[1] = 1;
    r_vector[2] = 2;
    r_vector[3] = r_vector[4] = 3;
    

    int64_t p = 23;
    int64_t precision = 6; 
    ZZ mod = power_ZZ(p, precision);
    zz_p::init(conv<long>(mod));

    map< Vec<int64_t>, zz_p, vi64less> f;
    Mat<int64_t> AP;
    Vec<int64_t> bP;
    stringstream buffer;
    buffer << input;
    buffer >> p;
    buffer >> f;
    buffer >> AP;
    buffer >> bP;
    
    const char charpoly_prec_str[] = "[47 44 42 40 38 36 34 32 30 28 26 24 22 20 18 16 14 12 10 8 6 4 3]";
    buffer << charpoly_prec_str;
    buffer >> charpoly_prec;

    int64_t verbose = 10;

    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    dr<zz_p> aa(p, f, AP, bP, verbose, false);
    int64_t weight = aa.n - 1;
    zeta_and_frob_Fp(zeta, F, hodge_numbers, r_vector, input, 100);
    Mat<zz_p> Fp;
    aa.frob_matrix(Fp, N);
    F = conv< Mat<ZZ> >(Fp);
    print(F);
    charpoly_frob(zeta, F, charpoly_prec, weight, p, 1);
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    return 0;
}

