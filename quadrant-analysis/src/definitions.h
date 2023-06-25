/************************************************************************
 * Quadrant analysis in base turbulent channel flow                     *
 * Amirreza Rastegari                                                   *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include "mpi.h"

/***************************************************************************************
 * Adjust these variables for the computations                                         *
****************************************************************************************/

#define zdv  221                               // channel full height (2h) in lattice units
#define ydv  256                               // channel width in lattice units
#define xdv  512                               // channel length in lattice units
#define CFL 0.1                                // CFL number (used in unit conversions)
#define Rebs 7200.                             // bulk Reynolds number Reb=(Ub x 2h)/nu
#define ub 0.0866034                           // nominal bulk velocity in lattice units (used in unit conversions)

// This is the factor that defines the vorticity, either 1. or 0.5
#define fac 1.

#define TSTR 600                               // beginning time of the computations
#define TEND 900                               // end time of the computations
#define DT 1                                   // Time step (i.e. calculate every DT time)

#define SLPTXT "base"
#define TEXT "Base Flow"

/***************************************************************************************
 * Don't touch from here on                                                            *
****************************************************************************************/

#define Fs(a,i,j,k,xd,yd,zd) (*(a+(i*yd*zd)+(j*zd)+k))
#define Fb(a,i,j,k,b,xd,yd,zd,d) (*(a+((i+2)*yd*zd*d)+(j*zd*d)+(k*d)+b))
#define Ts(i,j,k,xd,yd,zd) ((i*yd*zd)+(j*zd)+k)

#define cs 7
#define bs 12

#define MPDP MPI_DOUBLE

typedef double DP;
#define PI 3.1415

typedef struct pointer_to_arrays
{
	int * s;
	DP * vel;
	DP * den;
} POINTER;

typedef struct parallel_data
{
	int myrank;
	int numproc;
	int right;
	int left;
} PDATA;

extern DP w[19];
extern int e0[19];
extern int e1[19];
extern int e2[19];
extern DP E0[19];
extern DP E1[19];
extern DP E2[19];
extern int opposite[19];