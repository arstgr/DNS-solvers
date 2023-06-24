/* ****************************************************************************
 * by Amirreza Rastegari                                                      *
 * arstgri@gmail.com                                                         *                                               *
 *                                                                            *
 * To be used with postproc-spectra.c                                         *
 *                                                                            *
 *                                                                            *
 * Computes 1-dimensional energy spectra from DNS of turbulent channel flow   *
 *                                                                            *
 * Inputs:                                                                    *
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
 *    Euux.time, Evvx.time, Ewwx.time                                         *
 *    Euuy.time, Evvy.time, Ewwy.time                                         *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include "mpi.h"

#define xdv 512
#define ydv 256
#define zdv 221

#define TSTR 0
#define TEND 1
#define DT 1

#define Fb(a,i,j,k,b,xd,yd,zd,d) (*(a+(i*yd*zd*d)+(j*zd*d)+(k*d)+b))

typedef double DP;
#define PI 3.1415

typedef struct pointer_to_arrays
{
	DP * vel;
} POINTER;

typedef struct parallel_data
{
	int myrank;
	int numproc;
	int left;
	int right;
} PDATA;
