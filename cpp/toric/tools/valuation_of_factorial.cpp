#include "tools.h"

int64_t valuation_of_factorial(const int64_t n, const int64_t p)
{
    int64_t sum, tmp;
    sum = 0;
    tmp = n;
    while( tmp != 0 )
    {
        sum += tmp%p;
        tmp = tmp/p;
    }
    return (n - sum)/(p - 1);
}
