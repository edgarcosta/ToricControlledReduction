
// Copyright 2018 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"
#define verbose 0

using namespace std;

int main(int argc, char* argv[])
{
    if( argc != 4)
        printf("Usage: %s psi_num psi_den bound\n Note: skips p < 11\n", argv[0]);

    Vec<ZZ> zeta;
    Mat<ZZ> F;
    Vec<int64_t> hodge_numbers, r_vector, N, charpoly_prec;
    map< Vec<int64_t>, ZZ, vi64less> f;
    Mat<int64_t> AP;
    Vec<int64_t> bP;

    char input[] = "[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n";
    char keys_str[] = "[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]";

    stringstream buffer;
    buffer  << input;

    buffer >> f;
    buffer >> AP;
    buffer >> bP;
    Mat<int64_t> keys;
    buffer << keys_str;
    buffer >> keys;
    int64_t psi_num, psi_den, bound;
    psi_num = atol(argv[1]);
    psi_den = atol(argv[2]);
    bound = atol(argv[3]);

    //print(f);
    //print(bound);



    PrimeSeq s;
    s.reset(9);
    int64_t p = s.next(); //11
    printf("C2F2_%lld_%lld = \n", (long long int)psi_num, (long long int)psi_den);
    while (p <= bound) {
        // check for Tame ramification 
        int64_t C = (746496 % p); // 2^10 * 3^6;
        if (psi_num % p == 0 or psi_den %p == 0 or MulMod(C, PowerMod(psi_num, 12, p), p) == PowerMod(psi_den, 12, p)) {
            cout << "# " << p << " is tamely ramified" <<endl;
        }
        else {
            for(int64_t i = 1; i < keys.NumRows(); ++i) {
                f[keys[i]] = psi_den;
            }
            f[keys[0]] = - 12 * psi_num;
            zeta_and_frob_Fp(zeta, F, hodge_numbers, r_vector, p, f, AP, bP, verbose);
            cout << "[ "<< p << ", ";
            cout <<= zeta;
            cout << "],"<< endl;
        }
        p = s.next();
    }
    cout << "]" << endl;

    return 0;
}

