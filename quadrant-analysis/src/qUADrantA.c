#include "definitions.h"

POINTER qUADrantA(int xd, int yd, int zd, POINTER V, PDATA pd, long ts)
{
	int i,j,k,a,q;
	DP *UAVG,*VAVG,*WAVG,*DAVG,*DUMMY;
	DP *UTIN,*VTIN,*WTIN,*DTIN,*uw;
	double u,v,w;
	int extx,exty,num;
	char fn[50];
	MPI_Status status;
	FILE *sv;
	
	double quadrant[4][zd],dquadrant[4][zd];;
	int nquadrant[4][zd];
	double *DUMMYq;
	int *iDUMMYq;

	extx=xd-1;
	exty=yd;
	num=exty*extx;

	UAVG=(DP *)calloc(zd,sizeof(DP));
	VAVG=(DP *)calloc(zd,sizeof(DP));
	WAVG=(DP *)calloc(zd,sizeof(DP));
	DAVG=(DP *)calloc(zd,sizeof(DP));
	DUMMY=(DP *)calloc(zd,sizeof(DP));

	UTIN=(DP *)calloc(zd,sizeof(DP));
	VTIN=(DP *)calloc(zd,sizeof(DP));
	WTIN=(DP *)calloc(zd,sizeof(DP));
	DTIN=(DP *)calloc(zd,sizeof(DP));

	uw=(DP *)calloc(zd,sizeof(DP));
	
	DUMMYq=(DP *)calloc(4*zd,sizeof(DP));
	iDUMMYq=(int *)calloc(4*zd,sizeof(DP));

	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=1;k<(zd-1);k++)
			{
					*(UAVG+k) += Fb(V.vel,i,j,k,0,xd,yd,zd,3);
					*(VAVG+k) += Fb(V.vel,i,j,k,1,xd,yd,zd,3);
					*(WAVG+k) += Fb(V.vel,i,j,k,2,xd,yd,zd,3);
					*(DAVG+k) += (1./3.)*Fb(V.den,i,j,k,0,xd,yd,zd,1);
			}
		}
	}
	for (k=0;k<zd;k++)
	{
		*(UAVG+k) /= (double)(extx*exty*pd.numproc);
		*(VAVG+k) /= (double)(extx*exty*pd.numproc);
		*(WAVG+k) /= (double)(extx*exty*pd.numproc);
		*(DAVG+k) /= (double)(extx*exty*pd.numproc);
	}

	MPI_Reduce(UAVG,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(UAVG,DUMMY,zd*sizeof(DP));

	MPI_Reduce(VAVG,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(VAVG,DUMMY,zd*sizeof(DP));

	MPI_Reduce(WAVG,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(WAVG,DUMMY,zd*sizeof(DP));

	MPI_Reduce(DAVG,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(DAVG,DUMMY,zd*sizeof(DP));
	
	for (i=-1;i<(extx+1);i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=1;k<(zd-1);k++)
			{
				Fb(V.vel,i,j,k,0,xd,yd,zd,3) -= (*(UAVG+k));
				Fb(V.vel,i,j,k,1,xd,yd,zd,3) -= (*(VAVG+k));
				Fb(V.vel,i,j,k,2,xd,yd,zd,3) -= (*(WAVG+k));
			}
		}
	}
	
	for (a=0;a<4;a++)
	{
		for (k=0;k<zd;k++)
		{
			quadrant[a][k] = 0.;
			dquadrant[a][k] = 0.;
			nquadrant[a][k] = 0;
		}
	}

	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=1;k<(zd-1);k++)
			{
				u=Fb(V.vel,i,j,k,0,xd,yd,zd,3);
				v=Fb(V.vel,i,j,k,1,xd,yd,zd,3);
				w=Fb(V.vel,i,j,k,2,xd,yd,zd,3);
				
				*(uw+k) += u*w;

				if (u>0. && w>0.)
				{
				  quadrant[0][k] += u*w;
				  nquadrant[0][k] ++;
				}
				else if (u<0. && w>0.)
				{
				  quadrant[1][k] += u*w;
				  nquadrant[1][k] ++;
				}
				else if (u<0. && w<0.)
				{
				  quadrant[2][k] += u*w;
				  nquadrant[2][k] ++;
				}
				else if (u>0. && w<0.)
				{
				  quadrant[3][k] += u*w;
				  nquadrant[3][k] ++;
				}
			}
		}
	}
	
	MPI_Reduce(uw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uw,DUMMY,zd*sizeof(DP));
	
	MPI_Allreduce(&quadrant[0][0],DUMMYq,4*zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(&quadrant[0][0],DUMMYq,4*zd*sizeof(DP));
	
	MPI_Allreduce(&nquadrant[0][0],iDUMMYq,4*zd,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
	memcpy(&nquadrant[0][0],iDUMMYq,4*zd*sizeof(int));
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (k=1;k<(zd-1);k++)
	{
		dquadrant[0][k] = quadrant[0][k]/((double)(extx*exty*pd.numproc));
		dquadrant[1][k] = quadrant[1][k]/((double)(extx*exty*pd.numproc));
		dquadrant[2][k] = quadrant[2][k]/((double)(extx*exty*pd.numproc));
		dquadrant[3][k] = quadrant[3][k]/((double)(extx*exty*pd.numproc));
	}
	
	for (k=1;k<(zd-1);k++)
	{
		*(uw+k) /= (double)(extx*exty*pd.numproc);
		quadrant[0][k] /= (double)(nquadrant[0][k]);
		quadrant[1][k] /= (double)(nquadrant[1][k]);
		quadrant[2][k] /= (double)(nquadrant[2][k]);
		quadrant[3][k] /= (double)(nquadrant[3][k]);
	}

	if (!pd.myrank)
	{	
		sprintf(fn,"quadrant.%.3d.%.4d",pd.myrank,(int)ts);
		sv=fopen(fn,"wb");
		fwrite(&quadrant[0][0],sizeof(DP),4*zd,sv);
		fwrite(uw,sizeof(DP),zd,sv);
		fclose(sv);
		
		sprintf(fn,"dquadrant.%.3d.%.4d",pd.myrank,(int)ts);
		sv=fopen(fn,"wb");
		fwrite(&dquadrant[0][0],sizeof(DP),4*zd,sv);
		fwrite(uw,sizeof(DP),zd,sv);
		fclose(sv);
	}

	free(UAVG);
	free(VAVG);
	free(WAVG);
	free(DAVG);
	free(DUMMY);
	
	free(UTIN);
	free(VTIN);
	free(WTIN);
	free(DTIN);

	free(uw);
	
	free(iDUMMYq);
	free(DUMMYq);

return V;
}	