/***************************************************************************
 *   Amirreza Rastegari                                                    *
 *   arstgr@gmail.com                                                      *
 *                                                                         *
 *   Parallel LBM-BGK code with constant flow rate Q for DNS of turbulent  *
 *   channel flows                                                         *
 *   This program uses a new formulation for the constant flow rate        *
 *   implementation                                                        *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include "mpi.h"

/**************************************************************************
 *   Input Values of the program are listed here                          *
 *   zdv : Number of points in Z direction (Wall normall direction)       *
 *         a +2 is added to this value to conside the solid walls         *
 *   ydv : Number of points in Y direction (streamwise direction)         *
 *   xdv : Number of points in X direction (streamwise direction)         *
 *   CFL : The value of the CFL number for the program                    *
 *   TSTR: Time at which the computations start                           *
 *   TEND: Time at which computations finish                              *
 **************************************************************************/

#define Rebs 100.000
#define ub 0.1
#define zdv  25
#define ydv  8
#define xdv  16
#define TSTR 420
#define TEND 440
#define CFL  0.1

/**************************************************************************
 * Defining the slip pattern on the walls                                 *
 * PSLIP: slip with posts                                                 *
 * RSLIP: slip with ridges                                                *
 **************************************************************************/
//#define PSLIP
//#define RSLIP

/**************************************************************************
 * Defining the number of iterations                                      *
 * This value must be used whenever a predetermined number of iterations  *
 * is desired, otherwise it must be commented out                         *
 **************************************************************************/
//#define ITER 3

/**************************************************************************
 * Uses a 2D cartesian topology of a x b processors (x,y)                 *
 * Note: a x b = Number of Cores                                          *
 **************************************************************************/
#define XDIM 1
#define YDIM 1

/**************************************************************************
 * The input parameters of Iterpolation                                   *
 * spxe is the number of points on streamwise direction of initial        *
 * velocity field                                                         *
 * spye is the number of points in spanwise direction and spze is the     *
 * number of points in wall normal direction                              *
 * PARMT is the dimensional coefficient for bringing the velocity fields  *
 * it can be defined as u_tau_LB/u_tau_SP * V   where V is the center     *
 * line velocity from which u_tau was determined                          *
 **************************************************************************/

#define spxe 512
#define spye 280
#define spze 128
#define PARMT 0.1
#define FACTOR 1.

/**************************************************************************
 * Cahnnel Sizes                                                          *
 * Lx/h = 2xpi/alpha is the channel size in the streamwise direction      *
 * Ly/h = 2xpi/beta is the size in the spanwise direction                 *
 * Lz/h=2 is the size in wall normal direction                            *
 **************************************************************************/

#define ALPHA 0.339
#define BETA 0.678

/**************************************************************************
 * Code starts from here                                                  *
 **************************************************************************/

#define Fs(a,i,j,k,xd,yd,zd) (*(a+(i*yd*zd)+(j*zd)+k))
#define Fb(a,i,j,k,b,xd,yd,zd,d) (*(a+(i*yd*zd*d)+(j*zd*d)+(k*d)+b))
//#define Fb(a,i,j,k,b,xd,yd,zd,d) (*(a+(b*xd*yd*zd)+(i*yd*zd)+(j*zd)+k))
#define Ts(i,j,k,xd,yd,zd) ((i*yd*zd)+(j*zd)+k)

#define Vsp(a,i,j,k,b) (*(a+j*spx*(spz+1)*3+k*spx*3+i*3+b))
#define Vlb(a,i,j,k,b) (*(a+j*xdt*zd*3+i*zd*3+k*3+b))
#define Vspp(a,i,j,k) (*(a+j*spx*(spz+1)+k*spx+i))
#define Vlbp(a,i,j,k) (*(a+j*xdt*zd+i*zd+k))
#define VLB(a,i,j,k,b) (*(a+j*xd*zd*3+i*zd*3+k*3+b))
#define PLB(a,i,j,k) (*(a+j*xd*zd+i*zd+k))

#define N1(x,y,z) ((1./8.)*(1-x)*(1-y)*(1-z))
#define N2(x,y,z) ((1./8.)*(1+x)*(1-y)*(1-z))
#define N3(x,y,z) ((1./8.)*(1+x)*(1+y)*(1-z))
#define N4(x,y,z) ((1./8.)*(1-x)*(1+y)*(1-z))
#define N5(x,y,z) ((1./8.)*(1-x)*(1-y)*(1+z))
#define N6(x,y,z) ((1./8.)*(1+x)*(1-y)*(1+z))
#define N7(x,y,z) ((1./8.)*(1+x)*(1+y)*(1+z))
#define N8(x,y,z) ((1./8.)*(1-x)*(1+y)*(1+z))

#define STRIPE_COUNT "128"       /* must be an ascii string */
#define STRIPE_SIZE "1048576"    /* must be an ascii string */
#define CB_NODES_FCT "64"       /* max number of cores used in collective write operations */

#define cs 7
#define bs 12

#define MPDP MPI_DOUBLE

typedef double DP;
DP PI;

typedef struct pointer_to_arrays
{
	int * s;
	DP * f;
	DP * ftemp;
	DP force;
	DP * rightbufs;
	DP * rightbufr;
	DP * leftbufs;
	DP * leftbufr;
	
	DP * upbufs;
	DP * upbufr;
	DP * dwbufs;
	DP * dwbufr;
	
	DP * urbufs;
	DP * dlbufr;
	DP * dlbufs;
	DP * urbufr;
	
	DP rho;
	MPI_Request req[12];
} POINTER;

typedef struct parallel_data
{
	int myrank;
	int numproc;
	int right;
	int left;
	int up;
	int dw;
	int ur;
	int dl;
	int ul;
	int dr;
} PDATA;

extern const DP w[19];
extern const DP BI[19];

extern int e0[19];
extern int e1[19];
extern int e2[19];

extern DP E0[19];
extern DP E1[19];
extern DP E2[19];
extern int opposite[19];
extern int mirror[19];

#ifndef isnan
__inline int isnan(double var)
{
	return var!=var;
}
#endif