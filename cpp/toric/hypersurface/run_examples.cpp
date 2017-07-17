// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

bool run_examples(char examples[][3][buffer_length], const int64_t &examples_length, const int64_t &verbose )
{
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
    bool val = true;
    for(int64_t i = 0; i < examples_length; i++)
    {
        user_time = get_cpu_time();
        get_timestamp(&wtime1);
        cout << "Testing: "<< examples[i][0] <<endl;
        if( test_Fp( examples[i][1],  examples[i][2], verbose))
            cout <<"PASS"<<endl;
        else
        {
            cout <<"FAIL"<<endl;
            val = false;
        }
        get_timestamp(&wtime2);
        wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
        user_time = get_cpu_time() - user_time;
        printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );
    }
    return val;
}
