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
#include "definitions.h"
#include "fun_defs.h"

int main(int argc, char *argv[])
{
	
	int xd,yd,zd,TM;
	int i,j,k,q,a,b,c,nut;
	int cntr;
	long ss,up;
	POINTER V;
	PDATA pd;
	DP Fx=0.,Fy=0.,Fz=0.;
	DP Gx,Gy,Gz;
	FILE *mysave,*tm;
	FILE *sv;
	char fn[60],fn2[60];
	char *ch;
	MPI_Status status;
	int tst,ttt;
	DP U,Q,ux,dens,mdt;
	DP wup,wdo,dummy,davg,mu;
	DP pgrad,dl,dr,tau;
	DP Reb, dplus,Lplus,ubulk;
	int slp,extx,exty;

	const double ut=0.0400257;
	const double Ret=216.13878;
	Reb=Rebs;
	ubulk=ub;
	
	zd=zdv+2;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &(pd.numproc));
	MPI_Comm_rank(MPI_COMM_WORLD, &(pd.myrank));
	
	if (pd.myrank != 0)
		pd.left = pd.myrank - 1;
	else
		pd.left = pd.numproc - 1;
	
	if (pd.myrank != (pd.numproc-1))
		pd.right = pd.myrank + 1;
	else
		pd.right = 0;

	yd=ydv;
	xd=(xdv/pd.numproc);

	Fx=CFL*CFL*ut*ut*2./((double)(zd-2));
	tau=(3.*CFL*ut*((double)(zd-2))/(2.*Ret))+0.5;
	tau=1./tau;
	Fy=0.;
	Fz=0.;
	Gx=Fx;
	Gy=0.;
	Gz=0.;

	V.vel=(DP *)calloc(xd*yd*zd*3,sizeof(DP));
	V.s=(int *)calloc(xd*yd*zd,sizeof(int));
	V.velrecv=(DP *)calloc(xd*yd*zd*3,sizeof(DP));

	V.ycorl = (DP *)calloc(yd*zd*3,sizeof(double));
	V.xcorl = (DP *)calloc(xdv*zd*3,sizeof(double));
	V.xcorlsl = (DP *)calloc(xdv*zd*3,sizeof(double));
	V.xcorlnsl = (DP *)calloc(xdv*zd*3,sizeof(double));

	if (V.vel==NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.vel: %s\n",strerror(errno));
	}
	if (V.s==NULL)
        {
                fprintf(stderr,"\nUnable to allocate memory for V.s: %s\n",strerror(errno));
        }
	if (V.velrecv==NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.velrecv: %s\n",strerror(errno));
	}
	if (V.ycorl==NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.ycorl: %s\n",strerror(errno));
	}
	if (V.xcorl==NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.xcorl: %s\n",strerror(errno));
	}

//	slip_ridges_init(xd, yd, zd, V, pd);

	for (tst=TSTR;tst<=TEND;tst+=DT)
	{
		V=reading(xd, yd, zd, V, pd, tst, Fx, ubulk);
		V=mean_subtract(xd, yd, zd, V, pd, tst, Fx);
		MPI_Barrier(MPI_COMM_WORLD);
		V=ycorrel(xd,yd,zd,V,pd,tst,Fx);
		V=xcorrel(xd,yd,zd,V,pd,tst,Fx);
//		V=xcorrel_slp_nslp(xd,yd,zd,V,pd,tst,Fx);
		
		if (pd.myrank==0)
		{
			sprintf(fn,"xcorl.%.4d",tst);
			sv=fopen(fn,"wb");
			fwrite(V.xcorl,sizeof(double),(xdv*zd*3),sv);
			fclose(sv);
/*
			sprintf(fn,"xcorlsl.%.4d",tst);
			sv=fopen(fn,"wb");
			fwrite(V.xcorlsl,sizeof(double),(xdv*zd*3),sv);
			fclose(sv);

			sprintf(fn,"xcorlnsl.%.4d",tst);
			sv=fopen(fn,"wb");
			fwrite(V.xcorlnsl,sizeof(double),(xdv*zd*3),sv);
			fclose(sv);
*/
			sprintf(fn,"ycorl.%.4d",tst);
			sv=fopen(fn,"wb");
			fwrite(V.ycorl,sizeof(double),(yd*zd*3),sv);
			fclose(sv);
		}
		MPI_Barrier(MPI_COMM_WORLD); 
	}


	free(V.vel);
	free(V.velrecv);
	free(V.xcorl);
	free(V.xcorlsl);
	free(V.xcorlnsl);
	free(V.ycorl);

	MPI_Finalize();

	return(0);
}


