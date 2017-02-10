// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "../solve_system.h"
#include <iostream>
#include <sstream>

#define NTESTS 10
#define ZZbits 10

template<class R>
void test_input_plain(const char* input)
{
    stringstream buffer;
    buffer << input;
    Mat<R> T, U1, U2;
    R D1, D2;
    Vec<int64_t> B1, B2, initB;
    buffer >> T;
    buffer >> initB;
    buffer >> B1;
    buffer >> U1;
    buffer >> D1;
    solve_system(B2, U2, D2, T, initB);
    if( !( B1 == B2 and U1 == U2 and D1 == D2 ))
    {
        cout << "FAIL!\ninput\n:"; 
        cout << "T = \n"<< T << endl;
        cout << "initB = "<<initB<< endl;
        cout << "expected =?= output\n";
        cout << "B:" <<endl;
        cout << B1 << "=?=" << B2 << endl;
        cout << "D:" <<endl;
        cout << D1 << "=?=" << D2 << endl;
        cout << "U:" <<endl;
        cout << U1 << "\n=?=\n" << U2 << endl;
        assert(false);
    }
}

template<class R>
void test_input_padic(const char* input)
{
    stringstream buffer;
    buffer << input;
    Mat<R> T, U1, U2;
    R D1, D2;
    Vec<int64_t> B1, B2, initB;
    int64_t p, precision;
    ZZX f;
    buffer >> T;
    buffer >> initB;
    buffer >> p;
    buffer >> precision;
    buffer >> f;
    buffer >> B1;
    buffer >> U1;
    buffer >> D1;
    solve_system(B2, U2, D2, T, initB);
    if( !( B1 == B2 and U1 == U2 and D1 == D2 ))
    {
        cout << "FAIL!\ninput\n:"; 
        cout << "T = \n"<< T << endl;
        cout << "initB = "<<initB<< endl;
        cout << "expected =?= output\n";
        cout << "B:" <<endl;
        cout << B1 << "=?=" << B2 << endl;
        cout << "D:" <<endl;
        cout << D1 << "=?=" << D2 << endl;
        cout << "U:" <<endl;
        cout << U1 << "\n=?=\n" << U2 << endl;
        assert(false);
    }
}


// R = zz_p* or ZZ
template<class R>
void test_solve_system_right_inverse(const Mat<R> &T, const Vec<int64_t> &initB)
{
    Vec<int64_t> B;
    Mat<R> Unom;
    R Udenom;
    solve_system(B, Unom, Udenom, T, initB);
    int64_t i, j;
    int64_t ncols, nrows;
    ncols = T.NumCols();
    nrows = T.NumRows();

    // T + J, where J is the inclusion of B in range(T)
    Mat<R> S;
    S.SetDims(nrows, ncols + B.length());
    

    for(i = 0; i < nrows; i++)
        for(j = 0; j < ncols; j++)
            S[i][j] = T[i][j];
    for(i = 0; i < B.length(); i++)
        S[ B[i] ][ ncols + i] = R(1);

    Mat<R> SU;
    SU = S*Unom;
    
    
    Mat<R> DI;
    DI.SetDims(nrows, nrows);
    for(i = 0; i < nrows; i++)
        DI[i][i] = Udenom;
    
    if( DI != SU )
    {
        printf("U is not the right inverse!");
        cout << "T =\n" << T << endl;
        cout << "Unom =\n" << Unom << endl;
        cout << "Uden = " << Udenom << endl;
        assert(false);
    }
}

void random(ZZ &x)
{
    x = (int64_t)RandomBits_ulong(ZZbits) - (int64_t)RandomBits_ulong(ZZbits);
}


template<class R>
void test_solve_system_right_inverse_random(int64_t N)
{
    int64_t i, j, k, ncols, nrows;
    Vec<int64_t> initB;
    initB.SetLength(0);
    for(k = 0; k < N; k++)
    {
        ncols = RandomBnd(100);
        nrows = RandomBnd(100);
        Mat<R> T;
        T.SetDims(nrows, ncols);
        for(i = 0; i < nrows; i++)
            for(j = 0; j < ncols; j++)
                random(T[i][j]);

        test_solve_system_right_inverse(T, initB);
    }
}

int main()
{
    //test over ZZ
    {
        // DO OVER Fp once a and double check over ZZ with image
        Mat<ZZ> T, U1, U2;
        ZZ D1, D2;
        Vec<int64_t> B1, B2, initB;
        stringstream buffer;


        buffer << "[[-3  6 -1  1 -7][ 1 -2  2  3 -1][ 2 -4  5  8 -4]]\n[0 2]";
        buffer >> T;
        buffer >> B1;
        pivot_columns(B2, T);
        assert(B1 == B2);

        buffer << "[0 1]";
        buffer >> B1;
        pivot_columns(B2, transpose(T));
        assert(B1 == B2);

        buffer << "[[ 1 -3 -1  0 -1  5  0  0 -5  0][ 3  1 -1 -1 -1  1 -2 -1  2  0][ 1  1  1 -8  7 -1 -5  0  0 -1][ 0  0  0  0  0  1  1 -5 -3  1][ 1  0  1 -2 -1  0 23  1 -1  2]]\n[0 1 2 3 5]";
        buffer >> T;
        buffer >> B1;
        pivot_columns(B2, T);
        assert(B1 == B2);


        buffer << "[[  0   0   0   0   0   0   0   0   0  -1][  1  -2   1   1   2   2  -1 467  -3  -1][  0  -2   4 -31   1   1  -4   1   1   2][ -1   2  -1  -1   0   3  -7   3   4  -1]]\n[0 1 4 9]";
        buffer >> T;
        buffer >> B1;
        pivot_columns(B2, T);
        assert(B1 == B2);


        cout << "pivot_columns<ZZ> specific examples: PASS\n";



        const char input1[] = "[[1 3] [1 2] [2 5]]\n[]\n[2]\n[[ 2 -3  0] [-1  1  0] [ 1  1 -1]]\n-1";
        test_input_plain<ZZ>(input1);
        

        const char input2[] = "[[2 4 0] [1 2 0] [2 1 1]]\n[]\n[1]\n[[-1  0  4] [ 2  0 -2] [ 0  0  0] [-3  6  0]]\n6";
        test_input_plain<ZZ>(input2);
        
        const char input3[] = "[[  0   1   0  -1][  0  -2  -2   2][  0   1   4  -1][  0   1 -31  -1][  0   2   1   0][  0   2   1   3][  0  -1  -4  -7][  0 467   1   3][  0  -3   1   4][ -1  -1   2  -1]]\n[2 3 4]\n[2 3 4 6 7 8]\n[[  -26   -12     0     0     0    -4     0     0     0   -10][    8     1     0     0     0     2     0     0     0     0][  -10    -5     0     0     0     0     0     0     0     0][   -2     1     0     0     0     2     0     0     0     0][   30    20    10     0     0     0     0     0     0     0][ -320  -155     0    10     0     0     0     0     0     0][   -6     3     0     0    10    -4     0     0     0     0][  -46   -12     0     0     0    16    10     0     0     0][-3720  -465     0     0     0  -940     0    10     0     0][   42     4     0     0     0    -2     0     0    10     0]]\n10";
        test_input_plain<ZZ>(input3);
        cout <<"solve_system<ZZ> specific examples: PASS\n";

        test_solve_system_right_inverse_random<ZZ>(NTESTS);
        cout <<"solve_system<ZZ> random examples: PASS\n";
    }


    //test over Fp
    {
        zz_p::init(5UL);//
        const char input1[] = "[[3 1 0][1 2 0][2 1 1]]\n[]\n[1]\n[[1 0 4][3 0 3][0 0 0][3 1 0]]\n1";
        test_input_plain<zz_p>(input1);
        cout <<"solve_system<zz_p> specific examples: PASS\n";

        test_solve_system_right_inverse_random<zz_p>(NTESTS);
        cout <<"solve_system<zz_p> random examples: PASS\n";
    }
    //test over Zp
    {
        ZZ_p::init(ZZ(9765625UL));//5^10
        const char input1[] = "[[-2 -4  0][ 1  2  0][ 2  1  1]]\n[]\n5\n10\n[0 1]\n[1]\n[[8138021       0 3255209][3255208       0 3255208][      0       0       0][4882813       1       0]]\n1";
        test_input_padic<ZZ_p>(input1);
        cout <<"solve_system<ZZ_p> specific examples: PASS\n";


    }
    //test over Zq
    {
        ZZ_pX f;
        stringstream buffer;
        buffer << "[3 4 0 0 0 1]";
        buffer >> f;
        ZZ_pE::init(f);
        const char input1[] = "[[[9765623  2]  [9765621  4]  []]  [[2  9765624  0  0  9765624]  [4  9765623  0  0  9765623]  []]  [[4  9765623  0  0  9765623]  [2  9765624  0  0  9765624]  [2  9765624  0  0  9765624]]]\n[]\n5\n10\n[3 4 0 0 0 1]\n[1]\n[[[3865560  4679362  4679362  4679362  4679362]  []  [7678046  5413553  5413553  3998245  6828861]]  [[2034505  406901  406901  406901  406901]  []  [5926602  2176036  2176036  2883690  1468382]]  [[]  []  []]  [[1  4882813  4882813  4882813]  [1]  []]]\n[1]";
        test_input_padic<ZZ_pE>(input1);
        cout <<"solve_system<ZZ_pE> specific examples: PASS\n";
    }


    return 0;
}


/*
Some sage code to generate inputs:

p = 5
a = 1
precision = 10
padic = Padic(p, a, precision)
T = Matrix(padic.Zq, [[-2, -4, 0], [1, 2, 0], [2, 1, 1]])*(1 + padic.Zq.gen())
B, U, D = solve_system.solve_system_padic(padic, T)
I = '['+ T.str().replace('\n',"") +']'+ "\n" + Matrix(initB).str() + "\n" + str(padic.p) + "\n" + str(padic.prec) + "\n" + Matrix(padic.Fq.polynomial().list()).str()+ "\n"
O = Matrix(B).str() + "\n" + '['+ U.str().replace('\n',"") +']'+ "\n" + str(D)
I+O


def to_NTL(T):
    return str(([[ZZX(sum([ZZX(poly)*p**i for i,poly in enumerate(elt.list())])).list() for elt in row] for row in T.rows()])).replace(","," ")


p = 5
a = 5
precision = 10
padic = Padic(p, a, precision)
T = Matrix(padic.Zq, [[-2, -4, 0], [1, 2, 0], [2, 1, 1]])*(1 - padic.Zq.gen()) + Matrix(padic.Zq, [[0,0,0], [1, 2, 0], [2, 1, 1]])*(1 - padic.Zq.gen()**4)
B, U, D = solve_system.solve_system_padic(padic, T)
I = "";
I = to_NTL(T)+ "\n" + Matrix(initB).str() + "\n" + str(padic.p) + "\n" + str(padic.prec) + "\n" + Matrix(padic.Fq.polynomial().list()).str()+ "\n"
O = Matrix(B).str() + "\n" + to_NTL(U)+ "\n" + str(ZZ(D))
I+O
*/
