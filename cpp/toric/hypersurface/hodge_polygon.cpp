// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

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
 *     0          1         2        n-1    # slopes
 *  0 ----------- i ----------------- N -->  where N = dim PH^(n-1) (X)
 */
void hodge_polygon(Vec<int64_t> &hodge_polygon, Vec<int64_t> &slope, Vec<int64_t> hodge_numbers)
{
    int64_t dim_H = sum(hodge_numbers);
    hodge_polygon.SetLength( dim_H + 1, int64_t(0));
    slope.SetLength( dim_H + 1, int64_t(0));
    reverse(hodge_numbers.begin(), hodge_numbers.end());

    int64_t shift_vertical, shift_horizontal;
    int64_t i, k;
    shift_vertical = 0;
    shift_horizontal = 1;
    slope[0] = 0;
    hodge_polygon[0] = 0;
    int64_t s = 0;
    for(i = 0; i < hodge_numbers.length(); i++)
    {
        s += i * hodge_numbers[i];
        for(k = 0; k < hodge_numbers[i]; k++)
        {
            hodge_polygon[k + shift_horizontal] = i*(k + 1) + shift_vertical;
            slope[k + shift_horizontal] = i;
        }
        shift_horizontal += hodge_numbers[i];
        shift_vertical = hodge_polygon[shift_horizontal - 1];
    }
    assert_print(s, ==, shift_vertical);
}
