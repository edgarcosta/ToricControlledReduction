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
      "31 \n[ [0 2 1 0 0][1 1 1 0 0][0 0 0 0 0][2 0 0 1 0][0 1 0 0 2][0 0 0 0 3][1 0 2 0 0][1 0 0 2 0][0 0 0 3 0][1 2 0 0 0][0 1 2 0 0][0 3 0 0 0][0 0 2 0 0][0 0 3 0 0][2 0 1 0 0][2 1 0 0 0][0 2 0 0 1][3 0 0 0 0][0 0 1 0 0] ]\n[6 12 2 6 9 4 12 6 3 3 12 5 3 10 6 3 9 4 3] \n[ [-1 -1 -1 -1 -1][ 0  0  0  0  1][ 0  0  0  1  0][ 0  0  1  0  0][ 0  1  0  0  0][ 1  0  0  0  0] ] \n[3 0 0 0 0 0] \n",
      "[416787349596085853594397793593920583748646209303727133530079150721 -97932645670591822200019621870949081475631011551344027884614617 305721058284886021436065416870808787124758620867879379452491 -787745653377636049485922485371490854293020578025356906572 -126109926099037228765856477286719099382537513491612408 -459297337509500833174295182625928041455651714069348 358452656745187955130823503610245609876939423051 621666071358286429293831952150963596734199485 0 750078907107322975019655245951946959289 -1240825321931055376376600903146314242 478986647048005667091769674941513 -1343581057638164564072285203202 879454554192346983968633529 0 854614120410576490685 533578014824501451 -740309420498788 -220100912888 -1488715852 625611 -217 1]"
    },

  };


  return not run_examples(examples, examples_length, verbose);
}
