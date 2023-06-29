#include "definitions.h"

const DP w[19]={1./3.,1./18.,1./18.,1./18.,1./18.,1./18.,1./18.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.};
const DP BI[19]={0.,1./6.,1./6.,1./6.,1./6.,1./6.,1./6.,1./12.,1./12.,1./12.,1./12.,1./12.,1./12.,1./12.,1./12.,1./12.,1./12.,1./12.,1./12.};

int e0[19]={0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1, 1,-1,-1, 0, 0, 0, 0};
int e1[19]={0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1, 1,-1,-1};
int e2[19]={0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1};

DP E0[19]={0.,1.,-1.,0.,0.,0.,0.,1.,-1.,-1.,1.,1.,1.,-1.,-1.,0.,0.,0.,0.};
DP E1[19]={0.,0.,0.,1.,-1.,0.,0.,1.,1.,-1.,-1.,0.,0.,0.,0.,1.,1.,-1.,-1.};
DP E2[19]={0.,0.,0.,0.,0.,1.,-1.,0.,0.,0.,0.,1.,-1.,-1.,1.,1.,-1.,-1.,1.};
int opposite[19]={0,2,1,4,3,6,5,9,10,7,8,13,14,11,12,17,18,15,16};
int mirror[19]={0,1,2,3,4,6,5,7,8,9,10,12,11,14,13,16,15,18,17};