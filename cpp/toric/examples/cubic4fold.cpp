// Copyright 2018 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"
//this will take a while, output might be useful in the future
#define verbose 10

int main()
{


  int examples_length = 1;

  char examples[1][3][buffer_length] = {
    // f = 4*x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + 5*x1^3 + 6*x0^2*x2 + 12*x0*x1*x2 + 6*x1^2*x2 + 12*x0*x2^2 + 12*x1*x2^2 + 10*x2^3 + 6*x0^2*x3 + 6*x0*x3^2 + 3*x3^3 + 9*x1^2*x4 + 9*x1*x4^2 + 4*x4^3 + 3*x2^2 + 3*x2 + 2
    {
      "cubic 4 fold, cut out by x0**3 + x1**3 + x2**3 + (x0 + x1 + 2*x2)**3 + x3**3 + x4**3 + x5**3 + 2*(x0+x3)**3 + 3*(x1+x4)**3 + (x2+x5)**3, with H^4 \\cap H^2,2 = ZZ",
      "19 \n[ [0 2 1 0 0][1 1 1 0 0][0 0 0 0 0][2 0 0 1 0][0 1 0 0 2][0 0 0 0 3][1 0 2 0 0][1 0 0 2 0][0 0 0 3 0][1 2 0 0 0][0 1 2 0 0][0 3 0 0 0][0 0 2 0 0][0 0 3 0 0][2 0 1 0 0][2 1 0 0 0][0 2 0 0 1][3 0 0 0 0][0 0 1 0 0] ]\n[6 12 2 6 9 4 12 6 3 3 12 5 3 10 6 3 9 4 3] \n[ [-1 -1 -1 -1 -1][ 0  0  0  0  1][ 0  0  0  1  0][ 0  0  1  0  0][ 0  1  0  0  0][ 1  0  0  0  0] ] \n[3 0 0 0 0 0] \n",
      "[184144368549628275143663229532787625188711914273876985521 -268471159862411831380176745200156910903501843233528190 -223106227032475206133110868587388014601248069169137 4120151930424288201904171165048716797806981886780 -8559872431629407621684566132372680327853840485 15807705321568619799971497936052964594374590 8757731480093418171729361737425465149238 -48519287978356887378002004085459640716 201603689660762689936296415867560834 -651535469817792626756267271963493 1031318511781230908992904269036 2856837982773492822695025676 7913678622641254356495916 -38362708004493615301573 91086762179248789794 -168211933849028236 232980517796438 3226876977790 -13408076085 49521980 -20577 -190 1]"
    },

  };


  return not run_examples(examples, examples_length, verbose);
}
