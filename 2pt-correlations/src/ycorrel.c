#include "definitions.h"

POINTER ycorrel(int xd,int yd,int zd,POINTER V,PDATA pd,long tst,DP Fx)
{
	int i,j,k,a,b,c;
	double *ycor;
	int jto,jjto;

	ycor = (double *)calloc(3*yd*zd,sizeof(double));
	
	for (k=0;k<zd;k++)
	{
		for (j=0;j<yd;j++)
		{
			for (jto=0;jto<yd;jto++)
			{
				jjto = jto -j;
				if (jjto < 0)
					jjto += yd;

				for (i=0;i<xd;i++)
				{
					Fy(ycor,jjto,k,0,yd,zd,3) += Fb(V.vel,i,j,k,0,xd,yd,zd,3)*Fb(V.vel,i,jto,k,0,xd,yd,zd,3);
					Fy(ycor,jjto,k,1,yd,zd,3) += Fb(V.vel,i,j,k,1,xd,yd,zd,3)*Fb(V.vel,i,jto,k,1,xd,yd,zd,3);
					Fy(ycor,jjto,k,2,yd,zd,3) += Fb(V.vel,i,j,k,2,xd,yd,zd,3)*Fb(V.vel,i,jto,k,2,xd,yd,zd,3);
				}
			}
		}
	}
/*
	for (k=0;k<zd;k++)
	{
		for (j=0;j<yd;j++)
		{
			Fy(ycor,j,k,0,yd,zd,3) /= ((double)(yd*xdv));
			Fy(ycor,j,k,1,yd,zd,3) /= ((double)(yd*xdv));
			Fy(ycor,j,k,2,yd,zd,3) /= ((double)(yd*xdv));
		}
	}
*/	
	for (j=0;j<(yd*zd*3);j++)
		*(V.ycorl+j) = 0.;
			
	MPI_Reduce(ycor,V.ycorl,(3*yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
			
	free(ycor);
return V;
}