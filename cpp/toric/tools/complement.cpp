// Copyright 2017 Edgar Costa
// See LICENSE file for license details.


#include "tools.h"

// assumes B \subset [0, 1, ..., n - 1]
void complement(Vec<int64_t> &res, const int64_t n, const Vec<int64_t> B)
{
    int64_t i, j, position;
    res.SetLength(n - B.length());
    i = 0;
    position = 0;
    for(j = 0; j < B.length(); j++)
    {
        while(i < B[j]  && i < n)
        {
            res[position] = i;
            i++;
            position++;
        }
        i++;
    }
    while(i < n)
    {
        res[position] = i;
        i++;
        position++;
    }
    assert_print(position, ==, n - B.length());
}
