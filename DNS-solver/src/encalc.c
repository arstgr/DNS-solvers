#include "definitions.h"

POINTER encalc(int xd, int yd , int zd, POINTER V, DP Gx, int TM, long cntr,DP ubulk, MPI_Comm cart_grid)
{
	int i,j,k,a,q;
	double K,kt;
	FILE *en,*enl,*fh;
	double *uavg,*vavg,*wavg,*dummy;
	double ufluc,vfluc,wfluc;
	double dens,ux,uy,uz,dt,Re,umax,dplus,U;
	double momentum=0.,tmomentum=0.,tdens=0.,denst=0.;
	char fn[30];
	int rank, numproc;
	
	MPI_Comm_size(cart_grid, &numproc);
	MPI_Comm_rank(cart_grid, &rank);

	umax=1.;//uts*(2.5*log(Re)+5.5);

	dt=1./((double)(0.5*(zd-2.)/CFL));
	
	uavg=(double *)calloc(zd,sizeof(double));
	vavg=(double *)calloc(zd,sizeof(double));
	wavg=(double *)calloc(zd,sizeof(double));
	dummy=(double *)calloc(zd,sizeof(double));

	K=0.;
	kt=0.;
	for (i=0;i<(xd-1);i++)
	{
		for (j=0;j<(yd-1);j++)
		{
			for (k=1;k<(zd-1);k++)
			{
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<19;a++)
				{
					dens+=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					ux+= E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					uy+= E1[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					uz+= E2[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
				}
				momentum += ux;
				ux /= dens;
				uy /= dens;
				uz /= dens;
				denst += dens;
				(*(uavg+k)) += ux;
				(*(vavg+k)) += uy;
				(*(wavg+k)) += uz;
			}
		}
	}
	for (k=0;k<zd;k++)
	{
		(*(uavg+k)) /= ((double)((xd-1)*(yd-1)*XDIM*YDIM));
		(*(vavg+k)) /= ((double)((xd-1)*(yd-1)*XDIM*YDIM));
		(*(wavg+k)) /= ((double)((xd-1)*(yd-1)*XDIM*YDIM));
	}

	momentum /= ((double)((xd-1)*(yd-1)*(zd-2)*XDIM*YDIM));
	denst /= ((double)((xd-1)*(yd-1)*(zd-2)*XDIM*YDIM));

	MPI_Allreduce(uavg,dummy,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(uavg,dummy,zd*sizeof(DP));

	MPI_Allreduce(vavg,dummy,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(vavg,dummy,zd*sizeof(DP));

	MPI_Allreduce(wavg,dummy,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(wavg,dummy,zd*sizeof(DP));

	MPI_Allreduce(&momentum,&tmomentum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	momentum=tmomentum;

	MPI_Allreduce(&denst,&tdens,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	denst=tdens;

	V.rho=denst;	

	K=0.;
	U=0.;

	for (i=0;i<(xd-1);i++)
	{
		for (j=0;j<(yd-1);j++)
		{
			for (k=1;k<(zd-1);k++)
			{
				ux=0.;uy=0.;uz=0.;dens=0.;
				for (a=0;a<19;a++)
				{
					dens+=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					ux+= E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					uy+= E1[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					uz+= E2[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
				}
				ux /= dens;
				uy /= dens;
				uz /= dens;

				ufluc = ux-(*(uavg+k));
				vfluc = uy-(*(vavg+k));
				wfluc = uz-(*(wavg+k));

				K += 0.5*(ufluc*ufluc+vfluc*vfluc+wfluc*wfluc);
			}
		}
	}

	K /= ((double)((xd-1)*(yd-1)*(zd-2)*XDIM*YDIM));

	MPI_Allreduce(&K,&kt,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	K=kt;

	for (k=0;k<zd;k++)
		U += (*(uavg+k));

	U /= ((double)(zd-2));

	V.force=denst*ubulk-momentum;
//	V.force=denst*ubulk-momentum;
//	V.force=ubulk-momentum/denst;

	if (rank==0)
	{
			enl=fopen("energy.txt","a");
			fprintf(enl,"%f %.14f\n",((double)TM+dt*cntr),K);
			fclose(enl);

			fh=fopen("flowrate.txt","a");
			fprintf(fh,"%f %.12f\n",((double)TM+dt*cntr),U*((yd-1)*YDIM*(zd-2)));
			fclose(fh);

			fh=fopen("massrate.txt","a");
			fprintf(fh,"%f %.12f\n",((double)TM+(double)(dt*cntr)),(momentum*((double)((yd-1)*YDIM*(zd-2)))));
			fclose(fh);
	}


	free(dummy);
	free(uavg);
	free(vavg);
	free(wavg);

return V;
}