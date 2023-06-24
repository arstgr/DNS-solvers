#include "definitions.h"

POINTER xcorrel(int xd,int yd,int zd,POINTER V,PDATA pd,long tst,DP Fx)
{
	int i,j,k,a,b,c;
	double *xcor;
	int ito,iito;
	int num;
	MPI_Status status;
	int stag=100,rtag=100;

	memcpy(V.velrecv,V.vel,xd*yd*zd*3*sizeof(double));
	xcor = (double *)calloc(3*xdv*zd,sizeof(double));
	num=pd.myrank;
	
	for (a=0;a<pd.numproc;a++)
	{
		for (k=0;k<zd;k++)
		{
			for (j=0;j<yd;j++)
			{

				for (i=0;i<xd;i++)
				{
					for (ito=0;ito<xd;ito++)
					{
						iito = (num*xd + ito) - (pd.myrank*xd+i);
						if (iito<0)
							iito += xdv;
						 
						Fx(xcor,iito,k,0,xdv,zd,3) += Fb(V.vel,i,j,k,0,xd,yd,zd,3)*Fb(V.velrecv,ito,j,k,0,xd,yd,zd,3);
						Fx(xcor,iito,k,1,xdv,zd,3) += Fb(V.vel,i,j,k,1,xd,yd,zd,3)*Fb(V.velrecv,ito,j,k,1,xd,yd,zd,3);
						Fx(xcor,iito,k,2,xdv,zd,3) += Fb(V.vel,i,j,k,2,xd,yd,zd,3)*Fb(V.velrecv,ito,j,k,2,xd,yd,zd,3);
					}
				}
			}
		}
		MPI_Sendrecv_replace(V.velrecv, (3*xd*yd*zd), MPI_DOUBLE,pd.right, stag, pd.left, rtag, MPI_COMM_WORLD,&status);
		MPI_Sendrecv_replace(&num, 1, MPI_INT,pd.right, stag, pd.left, rtag, MPI_COMM_WORLD,&status);
	}
/*
	for (k=0;k<zd;k++)
	{
		for (i=0;i<xdv;i++)
		{
			Fx(xcor,i,k,0,xdv,zd,3) /= ((double)(yd*xdv));
			Fx(xcor,i,k,1,xdv,zd,3) /= ((double)(yd*xdv));
			Fx(xcor,i,k,2,xdv,zd,3) /= ((double)(yd*xdv));
		}
	}
*/			
	for (i=0;i<(xdv*zd*3);i++)
		*(V.xcorl+i) = 0.;

	MPI_Reduce(xcor,V.xcorl,(3*xdv*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	free(xcor);
return V;
}