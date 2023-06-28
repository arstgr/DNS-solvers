/************************************************************************
* Amirreza Rastegari                                                    *
* This should be used along with the budget-2d-revised                  *
* This code gives all of the statistics and budgets                     *
*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define Reb 7200.0000
#define ubulk 0.0866034
#define TSTR 500
#define TEND 800
// Time at which an instantaneous field is read and written
#define TINST 700
#define DT 1
#define TEXT "g<sup>+</sup>=56,w<sup>+</sup>=8"
// SLPF: 0 both, -1 no-slip, -2 slip
#define SLPF -3
#define SLPTXT "2DAVG"

// u_tau of the base channel flow with no-slip walls
#define utbase 0.005369 

//These are not in plus units, gplus = g+/delta+ !!! don't forget this
// width of the slip stripe
#define gplus 60
//width of the no-slip stripe
#define wplus 4
// size of the starting no-slip gap: half width of the no-slip stripe
#define wst 2
// width of the shear free area
#define shfw 64

typedef double DP;

#define xdv 2048
#define ydv 1024
#define zdv 221