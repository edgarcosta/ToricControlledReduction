
// Copyright 2018 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"
#define verbose 0

int main()
{

    int examples_length = 15;

    char examples[15][3][buffer_length] = {
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 251 = 2^8 - 5",
            "251 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[250058907189001 -353252214089 -3716113985 10206162 -58985 -89 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 509 = 2^9 - 3",
            "509 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[17390284781428441 111625490064943 363835479811 830613686 1404331 1663 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 1021 = 2^10 - 3",
            "1021 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[1132803161805372121 1141017400405050 194772803763 304392772 186843 1050 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 2039 = 2^11 - 9",
            "2039 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[-71862670932669131761 68448524227146360 -20675854993041 0 4973121 -3960 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 4093 = 2^12 - 3",
            "4093 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[-4701651857820438815449 741761249830962243 -190552118160103 0 11374447 -2643 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 8191 = 2^13 - 1",
            "8191 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[302010161517773079920641 -111720271586164924659 19248146778281775 -2464306827130 286889775 -24819 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 16381 = 2^14 - 3",
            "16381 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[-19321572190061729521304281 -2530321800382166363861 -149921787688268487 0 558706767 35141 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 32749 = 2^15 - 19",
            "32749 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[1233639479298456921244491001 23589323250194108972508 -648655336749212532 -28290325892378 -604808532 20508 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 65521 = 2^16 - 15",
            "65521 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[-79119421429263970001791209121 1191011361331792014508144 1328493693044639203 0 -309455683 -64624 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 131071 = 2^17 - 1",
            "131071 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[5070370291582725139136985169921 -6617604372820989664319382 28365273013510065867 -1579630508205868 1651101387 -22422 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 262139 = 2^18 - 5",
            "262139 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[-324481417228090561098992476241161 156477841628157718380152658 -1445644814364525391226 0 21037703306 -33138 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 524287 = 2^19 - 1",
            "524287 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[20768949750785132794698281184657409 -32349928100572518629675272311 -72811043957027987540593 245840516312049054 -264886045297 -428151 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 1048573 = 2^20 - 3",
            "1048573 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[-1329205178237316127846278475609539289 904772698614282742546403492979 -913044890036884128512599 0 830414241631 -748419 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 2097143 = 2^21 - 9",
            "2097143 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[85068401253498804942414792372842404849 26092794202314150547302646495389 -5232121010146346538750925 -6477448877457561486 -1189656795325 1348989 1]"
        },
        // f = x^3*y + y^4 + z^4 - 12*x*y*z + 1
        {
            "K3_C2F2, p = 4194301 = 2^22 - 3",
            "4194301 \n[ [ 1  0  0][ 3 -1 -1][ 0  0  0][ 0  0  2][ 0  2 -1] ]\n[-12 1 1 1 1] \n[ [ 1  0  0][ 1  3  0][-3 -3 -2][ 1  1  2] ] \n[0 0 4 0] \n",
            "[5444494505440932445170601091239952515801 1264895040620687336844356249944708 -24764922140050006711056828 -25983093852867318970 -1407724856028 4087108 1]"
        },

    };



    return not run_examples(examples, examples_length, verbose);
}