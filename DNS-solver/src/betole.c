#include "definitions.h"

void betole(void *pv, size_t n)
{
        char *p=pv,tmp;
        size_t lo,hi;
        for (lo=0,hi=n-1;hi>lo;lo++,hi--)
        {
                tmp=p[lo];
                p[lo]=p[hi];
                p[hi]=tmp;
        }
}