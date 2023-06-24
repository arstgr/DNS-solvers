/* ****************************************************************************
 * by Amirreza Rastegari                                                      *
 * arstgri@gmail.com                                                          *                                               *
 *                                                                            *
 * This code is used to calculate the 2pt velocity correlations in parallel   *
 * Formulation:							                                 	  *
 * Mathiu & Scott, An introduction to turbulent flows, Cambridge university   *
 * press, 2000, pp. 62						                             	  *
 * To be used with postproc-corl.c                                            *
 *                                                                            *
 *                                                                            *
 * Computes 2 point velocity auto correleations in DNS of turbulent channel   * 
 * flow                                                                       *
 *                                                                            *
 * Inputs:   vel.time  containing the instantaneous velocity field            *
 *    xd: dimension of velocity array in streamwise direction                 *
 *    yd: dimension of velocity array in spanwise direction                   *
 *    zd: dimension of velocity array in wall-normal direction (includes two  * 
 *        empty cells at the top and bottom in place for the walls)           *
 *    TSTR: starting time for the calculations                                *
 *    TEND: final time for the calculations                                   *
 *                                                                            *
 * Input files:                                                               *
 *    vel.time: velocity file with array indices running from                 *
 *    (u_x,u_y,u_z), then z index, then y index, then x index                 *
 *    example: vel.0001                                                       *
 *                                                                            *
 * Output files:                                                              *
 *    xcorl.time 2pt velocity autocorrelation in streamwise direction         *
 *    ycorl.time 2pt velocity autocorrelation in spanwise direction           *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include "mpi.h"


#define zdv  221
#define ydv  256
#define xdv  512
#define CFL 0.1
#define Rebs 7200.
#define ub 0.0866034

//Width of the slip stripe
#define gplus 0
//Width of the no-slip stripe
#define wplus 0
//half width of the no-slip stripe
#define gst 0

// This is the factor that defines the vorticity, either 1. or 0.5
#define fac 1.

#define TSTR 500
#define TEND 900
#define DT 1

#define Fs(a,i,j,k,xd,yd,zd) (*(a+(i*yd*zd)+(j*zd)+k))
#define Fb(a,i,j,k,b,xd,yd,zd,d) (*(a+(i*yd*zd*d)+(j*zd*d)+(k*d)+b))
#define Fy(a,j,k,b,yd,zd,d) (*(a+(k*yd*d)+(j*d)+b))
#define Fx(a,i,k,b,xd,zd,d) (*(a+(k*xd*d)+(i*d)+b))

#define MPDP MPI_DOUBLE

typedef double DP;
#define PI 3.1415

typedef struct pointer_to_arrays
{
	int * s;
	DP * vel;
	DP * velrecv;
	DP * xcorl;
	DP * xcorlsl;
	DP * xcorlnsl;
	DP * ycorl;
} POINTER;

typedef struct parallel_data
{
	int myrank;
	int numproc;
	int right;
	int left;
} PDATA;
