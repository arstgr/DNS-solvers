#include "definitions.h"

POINTER statistics(int xd, int yd, int zd, POINTER V, PDATA pd, long ts)
{
	int i,j,k;
	int a;
	DP u,v,w;
	DP UFLUC,VFLUC,WFLUC,DFLUC;
	DP *UAVG,*VAVG,*WAVG,*DUMMY,*DAVG;
	DP *UTIN,*VTIN,*WTIN,*REYS,*DTIN;
	DP dens, *zps;
	FILE *sv;
	int exty,extx;
	char fn[20];

	UAVG=(DP *)calloc(zd,sizeof(DP));
	VAVG=(DP *)calloc(zd,sizeof(DP));
	WAVG=(DP *)calloc(zd,sizeof(DP));
	DAVG=(DP *)calloc(zd,sizeof(DP));
	DUMMY=(DP *)calloc(zd,sizeof(DP));

	UTIN=(DP *)calloc(zd,sizeof(DP));
	VTIN=(DP *)calloc(zd,sizeof(DP));
	WTIN=(DP *)calloc(zd,sizeof(DP));
	REYS=(DP *)calloc(zd,sizeof(DP));
	DTIN=(DP *)calloc(zd,sizeof(DP));
	
	extx = xd-1;
	exty = yd-1;

	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=1;k<(zd-1);k++)
			{
					u=0.;v=0.;w=0.; dens=0.;
					for (a=0;a<19;a++)
					{
						dens += Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						u += E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						v += E1[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						w += E2[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					}

					*(UAVG+k) += (u/dens);
					*(VAVG+k) += (v/dens);
					*(WAVG+k) += (w/dens);
					*(DAVG+k) += dens;
			}
		}
	}

	for (k=0;k<zd;k++)
	{
		*(UAVG+k) /= (double)(extx*exty*XDIM*YDIM);
		*(VAVG+k) /= (double)(extx*exty*XDIM*YDIM);
		*(WAVG+k) /= (double)(extx*exty*XDIM*YDIM);
		*(DAVG+k) /= (double)(extx*exty*XDIM*YDIM);
	}
	
	MPI_Allreduce(UAVG,DUMMY,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(UAVG,DUMMY,zd*sizeof(DP));

	MPI_Allreduce(VAVG,DUMMY,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(VAVG,DUMMY,zd*sizeof(DP));

	MPI_Allreduce(WAVG,DUMMY,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(WAVG,DUMMY,zd*sizeof(DP));

	MPI_Allreduce(DAVG,DUMMY,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(DAVG,DUMMY,zd*sizeof(DP));
	
	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=1;k<(zd-1);k++)
			{
				u=0.;v=0.;w=0.; dens=0.;
				for (a=0;a<19;a++)
				{
					dens += Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					u += E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					v += E1[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					w += E2[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
				}

				UFLUC = (u/dens)-(*(UAVG+k));
				VFLUC = (v/dens)-(*(VAVG+k));
				WFLUC = (w/dens)-(*(WAVG+k));
				DFLUC = dens-(*(DAVG+k));
				*(UTIN+k) += UFLUC*UFLUC;
				*(VTIN+k) += VFLUC*VFLUC;
				*(WTIN+k) += WFLUC*WFLUC;
				*(REYS+k) += UFLUC*WFLUC;
				*(DTIN+k) += DFLUC*DFLUC;
			}
		}
	}

	for (k=0;k<zd;k++)
	{
		*(UTIN+k) /= (double)(extx*exty*XDIM*YDIM);
		*(VTIN+k) /= (double)(exty*extx*XDIM*YDIM);
		*(WTIN+k) /= (double)(exty*extx*XDIM*YDIM);
		*(REYS+k) /= (double)(extx*exty*XDIM*YDIM);
		*(DTIN+k) /= (double)(extx*exty*XDIM*YDIM);
	}

	MPI_Allreduce(UTIN,DUMMY,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(UTIN,DUMMY,zd*sizeof(DP));

	MPI_Allreduce(VTIN,DUMMY,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(VTIN,DUMMY,zd*sizeof(DP));

	MPI_Allreduce(WTIN,DUMMY,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(WTIN,DUMMY,zd*sizeof(DP));

	MPI_Allreduce(REYS,DUMMY,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(REYS,DUMMY,zd*sizeof(DP));

	MPI_Allreduce(DTIN,DUMMY,zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(DTIN,DUMMY,zd*sizeof(DP));

	for (k=0;k<zd;k++)
	{
		*(UTIN+k) = sqrt((*(UTIN+k)));
		*(VTIN+k) = sqrt((*(VTIN+k)));
		*(WTIN+k) = sqrt((*(WTIN+k)));
		*(DTIN+k) = sqrt((*(DTIN+k)));
	}


	if (pd.myrank==0)
	{
		sprintf(fn,"turb-fld.%.3d.%.3ld",pd.myrank,ts);
		sv=fopen(fn,"wb");
		fwrite(UAVG,sizeof(DP),zd,sv);
		fwrite(VAVG,sizeof(DP),zd,sv);
		fwrite(WAVG,sizeof(DP),zd,sv);
		fwrite(UTIN,sizeof(DP),zd,sv);
		fwrite(VTIN,sizeof(DP),zd,sv);
		fwrite(WTIN,sizeof(DP),zd,sv);
		fwrite(REYS,sizeof(DP),zd,sv);
		fwrite(DAVG,sizeof(DP),zd,sv);
		fwrite(DTIN,sizeof(DP),zd,sv);
		fclose(sv);
	
// Test print function
		zps = (double *)calloc(zd,sizeof(double));
		for (k=1;k<zd;k++)
		  *(zps+k) = k-0.5;
		*(zps+zd-1) = zd-2;
	
		for (k=0;k<zd;k++)
		  *(zps+k) /= ((double)(zd-2.));
		
		sv=fopen("uavg-temp.dat","w");
		for (k=0;k<zd;k++)
		  fprintf(sv,"%.12f %.12f\n",*(zps+k),*(UAVG+k));
		fclose(sv);
		
		free(zps);
// End of test print function
	}

	free(UAVG);
	free(VAVG);
	free(WAVG);
	free(DAVG);

	free(UTIN);
	free(VTIN);
	free(WTIN);
	free(REYS);
	free(DTIN);

	free(DUMMY);
	
//	fprintf(stderr,"in function statistics, c error handler is %s\n",strerror(errno));

	return V;
}