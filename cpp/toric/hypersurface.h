// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#ifndef HYPERSURFACE_H
#define HYPERSURFACE_H

#include "dr.h"

using namespace std;
using namespace NTL;

/*
 * Let X be a nondegenerate toric hypersurface defined over Fp of dimension `weight`
 * This is a minimal example on how to compute compute the interesting factor of the zeta function of X
 *
 * Input:
 *  - p, the characteristic of the field
 *  - f, the defining polynomial of X such that convex(supp(f)) = P
 *  - AP, the matrix associated to the half space representation of P,  w \in k*P <=>  AP * v + k *bP >= 0
 *  - bP,  the vector associated to the half space representation of P,  w \in k*P <=>  AP * v + k *bP >= 0
 *  - verbose, the verbose level, check dr.h for more details
 *
 *  Output:
 *  - zeta, the coefficients of the interesting factor of the zeta function of X, i.e., the characteristic polynomial of Frobenius acting on PH^(weight) (X)
 *  - F, a matrix representing a p-adic approximation of the p-power Frobenius acting on PH^(weight) (X)
 *  - hodge_numbers, the hodge numbers  [h_0, h_1,   ...,    h_n], the first h0 columns correspond to PH^(weight, 0)(X), the next h1 to PH^(weight-1, 1)(X), etc
 *  - r_vector, a column corresponding to basis element in PH^(weight - i, i)(X) is correct modulo p^(r_vector[i] + weight - i)
 */
void zeta_and_frob_Fp(Vec<ZZ> &zeta, Mat<ZZ> &F, Vec<int64_t> &hodge_numbers, Vec<int64_t>& r_vector, const int64_t &p, const map< Vec<int64_t>, ZZ, vi64less> &f, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const int64_t &verbose = 1);
/*
 * same as above, but the input is given through string in the following format:
 *      p
 *      f.keys()
 *      f.values()
 *      AP
 *      bP
 *
 *  for example:
 *  '17 \n[[0 2][0 0][3 0][2 0][0 1]]\n[1 12 16 1 1] \n[[ 0  1][ 1  0][-2 -3]] \n[0 0 6] \n'
 *  represents:
 *  y**2 + y - (x ** 3 - x**2 - 7820*x - 263580) over GF(17) aka http://www.lmfdb.org/EllipticCurve/Q/11/a/1
 */
void zeta_and_frob_Fp(Vec<ZZ> &zeta, Mat<ZZ> &F, Vec<int64_t> &hodge_numbers, Vec<int64_t>& r_vector, const char* input, const int64_t &verbose = 1);

/*
 * Given an input in format above, computes the zeta function
 * and checks if it maches the expected output
 */
bool test_Fp(const char *input, const char *output, const int64_t &verbose = 0);


/*
 * Same goal as above, but now over Fq. However, we have not implemented the Frobenius action on Zq. Thus some of the works is left to you.
 * Explicitly, we compute an approximation of the p-power frobenius, and you will need to deduce from it a q-power approximation and then lift the right characteristic polynomial, perhaps with charpoly_frob (if doing this over Sage, check https://github.com/edgarcosta/yellow/blob/master/charpoly_frobenius.py )
 *
 * Let X be a nondegenerate toric hypersurface defined over Fq of dimension `weight`
 *
 * Input:
 *  - p, the characteristic of the field
 *  - fE, the polynomial defining the extension Fq = Fp[x]/fE(X)
 *  - f, the defining polynomial of X such that convex(supp(f)) = P
 *  - frob, sigma(f)
 *  - AP, the matrix associated to the half space representation of P,  w \in k*P <=>  AP * v + k *bP >= 0
 *  - bP,  the vector associated to the half space representation of P,  w \in k*P <=>  AP * v + k *bP >= 0
 *  - verbose, the verbose level, check dr.h for more details
 *
 *  Output:
 *  - F, a matrix representing a p-adic approximation of the p-power Frobenius acting on PH^(weight) (X)
 *  - hodge_numbers, the hodge numbers  [h_0, h_1,   ...,    h_n], the first h0 columns correspond to PH^(weight, 0)(X), the next h1 to PH^(weight-1, 1)(X), etc
 *  - r_vector, a column corresponding to basis element in PH^(weight - i, i)(X) is correct modulo p^(r_vector[i] + weight - i)
 *  - charpoly_prec,  if we compute the charpoly of F F^\sigma .... F^(\sigma^(a - 1)), then the ith coefficient will be correct mod p^charpoly_prec[i]
 *
 *  Note:
 *   It's up to you to compute F F^\sigma .... F^(\sigma^(a - 1)) and then deduce the characteristic polynomial of Frob with charpoly_frob()
 */
void frob_Fq(Mat<ZZX> &F, Vec<int64_t> &hodge_numbers, Vec<int64_t>& r_vector, Vec<int64_t> &charpoly_prec, const int64_t p, const ZZX &fE, map< Vec<int64_t>, ZZX, vi64less> &f, map< Vec<int64_t>, ZZX, vi64less> &ffrob, const Mat<int64_t> &AP, const Vec<int64_t> &bP, const int64_t &verbose = 1);

/*
 * same as above, but the input is given through string in the following format:
 *      p
 *      fE
 *      f.keys()
 *      f.values()
 *      ffrob.keys()
 *      ffrob.values()
 *      AP
 *      bP
 *
 *  for example:
 *  '17 \n[1 1 1] \n[[0 2][0 0][3 0][2 0][0 1]]\n[[1] [12] [16] [1] [1] ] \n[[2 0][0 0][3 0][0 2][0 1]]\n[[1] [12] [16] [1] [1] ] \n[[ 0  1][ 1  0][-2 -3]] \n[0 0 6] \n
 *  represents:
 *  y**2 + y - (x ** 3 - x**2 - 7820*x - 263580) over GF(17^2) aka http://www.lmfdb.org/EllipticCurve/Q/11/a/1
 */
void frob_Fq(Mat<ZZX> &F, Vec<int64_t> &hodge_numbers, Vec<int64_t>& r_vector, Vec<int64_t> &charpoly_prec, const char* input, const int64_t &verbose = 1);

/*
 * Given an input in format above, computes the Frob approximation
 * and checks if it matches the expected output
 */
bool test_Fq(const char *input, const char *output, const int64_t &verbose = 0);


/*
 *  Computes the hodge polygon
 *
 *  Input:
 *   - hodge_numbers = [h_0, h_1,   ...,    h_n], see below
 *
 *  Output:
 *   - hodge_polygon, the hodge polygon, where hodge_polygon[i] corresponds to point in the graph below above i
 *      e.g hodge_polygon = [0 0 0 0 1 2 3 4 6 8 .... sum(i*h_{n-i})]
 *   - slope, slope[i] corresponds to slope at the point i, basically which h_{n-i} is making the hodge polygon grow
 *      e.g, slope = [0 0 0 0 1 1 1 1 ..... n-1 n-1 ]
 *
 *
 *  PH^(n-1) (X) = H^{n - 1, 0} + H^{n - 2, 1} + ... + H^{0, n - 1}
 *  write dim H^{n - 1 -i, i} = h_i
 *
 *                                     !
 *                                    !
 *
 *                                 .
 *                               .
 *                             .
 *
 *                         /
 *                       /
 *                     ~
 *                   -
 *                 -
 *               -
 *             -
 *           -
 *         ~
 *   -----
 *  h_{n-1} ~  h_{n-2} ~ h_{n-3}  ... h_0   # hodge numbers
 *     0          1         2       n -1    # slopes
 *  0 ----------- i ----------------- N -->  where N = dim PH^(n-1) (X)
 *
 */
void hodge_polygon(Vec<int64_t> &np, Vec<int64_t> &slope, Vec<int64_t> hodge_numbers);


/*
 * Input:
 *  - hodge_numbers, hodge_numbers[i] = dim PH^(weight - i, i)(X)
 *  - weight, weight of motive = dim(X) = n - 1
 *  - p, the characteristic of the field
 *  - a, q = p^a
 *
 * Output:
 *  - r_vector, the number of digits above the Hodge polygon necessary to deduce the correct characteristic polynomial, where r_vector[i] corresponds to PH^(weight - i, i)(X)
 *  - charpoly_prec, if compute the characteristic polynomial of a matrix representing an approximation of the Frobenius acting on PH^weight (X) with enough digits (according to r_vector), then the ith coefficient will be correct mod p^charpoly_prec[i]
 */
void relative_precision(Vec<int64_t> &r_vector, Vec<int64_t> &charpoly_prec, const Vec<int64_t> &hodge_numbers, const int64_t &weight, const int64_t &p, const int64_t &a);

/*
 * Input:
 *  - r_vector, corresponds to desired number digits above the above the Hodge polygon, where r_vector[i] corresponds to PH^(weight - i, i)(X)
 *  - weight, weight of motive = dim(X) = n - 1
 *  - p, the characteristic of the field
 *
 *  Output:
 *   - N, that determines how many terms of the Frobenius expansion we must use to assure the desired number of digits above the Hodge polygon
 */
void N_vector(Vec<int64_t> &N, const Vec<int64_t> &r_vector, const int64_t &weight, const int64_t &p);

// outputs =  max([ self.r_vector[m - 1] + max(factorial_p_adic(self.p*(m + self.N[m - 1] - 1) - 1, self.p)[0] - m + 1, 0) for m in range(1, self.n + 1) ])
int64_t working_precision(const Vec<int64_t> &N, const Vec<int64_t> &r_vector, const int64_t &p);

/*
 *  Input:
 *   - frob_matrix, a matrix over ZZ representing an approximation of Frob in H^(weight + 1) of the complement
 *   - charpoly_prec, if we compute the charpoly of frob_matrix, then the ith coefficient will be correct mod p^charpoly_prec[i]
 *   - weight, weight of motive = n - 1
 *   - p, the characteristic of the field
 *   - q, q = p^a and Frob is a lift of the qth power
 *
 *  Output:
 *   - a list of integers corresponding to the characteristic polynomial of the Frobenius action
 *
 */
void charpoly_frob(Vec<ZZ> &cp, const Mat<ZZ> &frob_matrix, const Vec<int64_t> &charpoly_prec, const int64_t &weight, const int64_t &p, const int64_t &a);

#define buffer_length 4096
bool run_examples(char examples[][3][buffer_length], const int64_t &examples_length, const int64_t &verbose = 0);

#endif
