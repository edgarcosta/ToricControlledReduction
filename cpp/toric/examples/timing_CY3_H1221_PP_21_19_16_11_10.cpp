
// Copyright 2018 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"
#define verbose 0

int main()
{

int examples_length = 15;

char examples[15][3][buffer_length] = {
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 251 = 2^8 - 5",
 "251 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 509 = 2^9 - 3",
 "509 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 1021 = 2^10 - 3",
 "1021 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 2039 = 2^11 - 9",
 "2039 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 4093 = 2^12 - 3",
 "4093 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 8191 = 2^13 - 1",
 "8191 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 16381 = 2^14 - 3",
 "16381 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 32749 = 2^15 - 19",
 "32749 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 65521 = 2^16 - 15",
 "65521 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 131071 = 2^17 - 1",
 "131071 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 262139 = 2^18 - 5",
 "262139 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 524287 = 2^19 - 1",
 "524287 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 1048573 = 2^20 - 3",
 "1048573 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 2097143 = 2^21 - 9",
 "2097143 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},
// f = y^7 + x^2*z*w + x*y*z*w + y^2*z*w + z^3*w + w^3 + x*z + y*z
{
"CY3_H1221_PP_21_19_16_11_10, p = 4194301 = 2^22 - 3",
 "4194301 \n[ [ 0  0 -3  2][ 2  0 -2  1][ 0  2 -2  1][ 0  1 -2  1][ 0  0  0  0][ 1  1 -2  1][ 0  7 -3  1][ 1  0 -2  1] ]\n[1 1 1 1 1 1 1 1] \n[ [  1   0   0   0][  1   1   2   3][  1   1   6   7][  1   1   6  11][  0   1   7  14][  0   1   7   7][  0   1   0   0][ -7  -6 -28 -42] ] \n[0 0 4 0 0 7 0 0] \n",
"[]"
},

};


    
    return not run_examples(examples, examples_length, verbose);
}
