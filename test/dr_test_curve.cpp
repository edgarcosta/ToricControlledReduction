// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "dr.h"
#include <iostream>
#include <sstream>
// tests if the hodge numbers match the expected hodge numbers for a sample of curves
// runs test_all() on a sample of curves

int main()
{
    {
        //elliptic curve 11.a
        Vec<int64_t> H;
        H.SetLength(3, 1);
        H[0] = 0;

        {
            const char ec[] = "17 \n[[0 2][0 0][3 0][2 0][0 1]]\n[1 12 16 1 1] \n[[ 0  1][ 1  0][-2 -3]] \n[0 0 6] \n";
            zz_pPush push(17);
            dr<zz_p> EC(ec, 0, false);
            assert_print(EC.cokernels_I_dimensions, ==, H);
            EC.test_all();
            cout << "Elliptic curve 11.a p = 17 over zz_p: PASS\n";
        }
        {
            const char ec[]  = "17 \n[[0 2][0 0][3 0][2 0][0 1]]\n[1 12 16 1 1] \n[[ 0  1][ 1  0][-2 -3]] \n[0 0 6] \n";
            ZZ_pPush push(ZZ(17));
            dr<ZZ_p> EC(ec, 0, false);
            assert_print(EC.cokernels_I_dimensions, ==, H);
            EC.test_all();
            cout << "Elliptic curve 11.a p = 17 over ZZ_p: PASS\n";

        }
        {
            const char ec[] = "0 \n[[0 2][0 0][3 0][2 0][0 1]]\n[1 12 16 1 1] \n[[ 0  1][ 1  0][-2 -3]] \n[0 0 6] \n";
            dr<ZZ> EC(ec, 0, false);
            assert_print(EC.cokernels_I_dimensions, ==, H);
            EC.test_all();
            cout << "Elliptic curve 11.a over ZZ: PASS\n";
        }
    }
    {
        //genus 2 curve 11664.a
        Vec<int64_t> H;
        H.SetLength(3, 2);
        H[0] = 0;
        {
            zz_pPush push(31);
            const char g2[] = "31 \n[[6 0][0 0][0 2]]\n[1 23 1] \n[[ 1  0][ 0  1][-1 -3]] \n[0 0 6] \n";
            dr<zz_p> G2(g2, 0, false);
            assert_print(G2.cokernels_I_dimensions, ==, H);
            G2.test_all();
            cout << "Genus 2 curve 11664.a p=31 over zz_p : PASS\n";
        }
        {
            ZZ_pPush push(ZZ(31));
            const char g2[] = "31 \n[[6 0][0 0][0 2]]\n[1 23 1] \n[[ 1  0][ 0  1][-1 -3]] \n[0 0 6] \n";
            dr<ZZ_p> G2(g2, 0, false);
            assert_print(G2.cokernels_I_dimensions, ==, H);
            G2.test_all();
            cout << "Genus 2 curve 11664.a p=31 over ZZ_p : PASS\n";

        }
        {
            const char g2[] = "0 \n[[6 0][0 0][0 2]]\n[1 23 1] \n[[ 1  0][ 0  1][-1 -3]] \n[0 0 6] \n";
            dr<ZZ> G2(g2, 0, false);
            assert_print(G2.cokernels_I_dimensions, ==, H);
            G2.test_all();
            cout << "Genus 2 curve 11664.a over ZZ : PASS\n";
        }
    }
    {
        //C_{3,4} curve
        Vec<int64_t> H;
        H.SetLength(3, 3);
        H[0] = 0;
        {
            zz_pPush push(43);
            const char c34[] = "43 \n[[0 3][0 0][4 0][1 0]]\n[42 1 1 1] \n[[ 0  1][ 1  0][-3 -4]] \n[0 0 12] \n";
            dr<zz_p> C34(c34, 0, false);
            assert_print(C34.cokernels_I_dimensions, ==, H);
            C34.test_all();
            cout << "C34 curve p=43 over zz_p : PASS\n";
        }
        {
            ZZ_pPush push(ZZ(43));
            const char c34[] = "43 \n[[0 3][0 0][4 0][1 0]]\n[42 1 1 1] \n[[ 0  1][ 1  0][-3 -4]] \n[0 0 12] \n";
            dr<ZZ_p> C34(c34, 0, false);
            assert_print(C34.cokernels_I_dimensions, ==, H);
            C34.test_all();
            cout << "C34 curve p=43 over ZZ_p : PASS\n";
        }
        {
            const char c34[] = "0 \n[[0 3][0 0][4 0][1 0]]\n[42 1 1 1] \n[[ 0  1][ 1  0][-3 -4]] \n[0 0 12] \n";
            dr<ZZ> C34(c34, 0, false);
            assert_print(C34.cokernels_I_dimensions, ==, H);
            C34.test_all();
            cout << "C34 curve over ZZ : PASS\n";
        }
    }
    {
        //plane curve, genus 3
        Vec<int64_t> H;
        H.SetLength(3, 3);
        H[0] = 0;
        {
            zz_pPush push(17);
            const char g3[] = "17 \n[[0 3][0 0][3 1][2 2][1 3][1 2][0 1][4 0][1 1][0 2][2 0][0 4]]\n[7 1 14 16 1 16 6 15 15 15 5 1] \n[[-1 -1][ 0  1][ 1  0]] \n[4 0 0] \n";
            dr<zz_p> G3(g3, 0, false);
            assert_print(G3.cokernels_I_dimensions, ==, H);
            G3.test_all();
            cout << "plane quartic curve p=17 over zz_p : PASS\n";
        }
        {
            ZZ_pPush push(ZZ(17));
            const char g3[] = "17 \n[[0 3][0 0][3 1][2 2][1 3][1 2][0 1][4 0][1 1][0 2][2 0][0 4]]\n[7 1 14 16 1 16 6 15 15 15 5 1] \n[[-1 -1][ 0  1][ 1  0]] \n[4 0 0] \n";

            dr<ZZ_p> G3(g3, 0, false);
            assert_print(G3.cokernels_I_dimensions, ==, H);
            G3.test_all();
            cout << "plane quartic curve p=17 over ZZ_p : PASS\n";
        }
        {
            const char g3[] = "0 \n[[0 3][0 0][3 1][2 2][1 3][1 2][0 1][4 0][1 1][0 2][2 0][0 4]]\n[7 1 14 16 1 16 6 15 15 15 5 1] \n[[-1 -1][ 0  1][ 1  0]] \n[4 0 0] \n";
            dr<ZZ> G3(g3, 0, false);
            assert_print(G3.cokernels_I_dimensions, ==, H);
            G3.test_all();
            cout << "plane quartic curve ZZ : PASS\n";
        }
    }

    {
        //genus 9 curve
        Vec<int64_t> H;
        H.SetLength(3, 9);
        H[0] = 0;
        {
            zz_pPush push(37);
            const char g9[] = "37 \n[[-1 -1][ 0  4][ 3  0]]\n[1 1 36] \n[[ 5 -1][-1  4][-4 -3]] \n[4 3 12] \n";
            dr<zz_p> G9(g9, 0, false);
            assert_print(G9.cokernels_I_dimensions, ==, H);
            G9.test_all();
            cout << "Genus 9 curve p=37 over zz_p : PASS\n";
        }
        {
            ZZ_pPush push(ZZ(37));
            const char g9[] = "37 \n[[-1 -1][ 0  4][ 3  0]]\n[1 1 36] \n[[ 5 -1][-1  4][-4 -3]] \n[4 3 12] \n";
            dr<ZZ_p> G9(g9, 0, false);
            assert_print(G9.cokernels_I_dimensions, ==, H);
            G9.test_all();
            cout << "Genus 9 curve p=37 over ZZ_p : PASS\n";
        }
        {
            const char g9[] = "0 \n[[-1 -1][ 0  4][ 3  0]]\n[1 1 36] \n[[ 5 -1][-1  4][-4 -3]] \n[4 3 12] \n";
            dr<ZZ> G9(g9, 0, false);
            assert_print(G9.cokernels_I_dimensions, ==, H);
            G9.test_all();
            cout << "Genus 9 curve ZZ : PASS\n";
        }
    }

    {
        //hyperelliptic curve of genus 10
        Vec<int64_t> H;
        H.SetLength(3, 10);
        H[0] = 0;
        {
            zz_pPush push(13);
            const char h10[] = "13 \n[[22  0][ 1  0][ 0  0][ 0  2]]\n[1 1 1 12] \n[[  0   1][  1   0][ -1 -11]] \n[0 0 22] \n";
            dr<zz_p> H10(h10, 0, false);
            assert_print(H10.cokernels_I_dimensions, ==, H);
            H10.test_all();
            cout << "Hyperelliptic curve of genus 10 p=13 over zz_p : PASS\n";
        }
        {
            ZZ_pPush push(ZZ(13));
            const char h10[] = "13 \n[[22  0][ 1  0][ 0  0][ 0  2]]\n[1 1 1 12] \n[[  0   1][  1   0][ -1 -11]] \n[0 0 22] \n";
            dr<ZZ_p> H10(h10, 0, false);
            assert_print(H10.cokernels_I_dimensions, ==, H);
            H10.test_all();
            cout << "Hyperelliptic curve of genus 10 p=13 over ZZ_p : PASS\n";
        }
        {
            const char h10[] = "0 \n[[22  0][ 1  0][ 0  0][ 0  2]]\n[1 1 1 12] \n[[  0   1][  1   0][ -1 -11]] \n[0 0 22] \n";
            dr<ZZ> H10(h10, 0, false);
            assert_print(H10.cokernels_I_dimensions, ==, H);
            H10.test_all();
            cout << "Hyperelliptic curve of genus 10 ZZ : PASS\n";
        }
    }



    return 0;
}
