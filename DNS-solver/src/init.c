#include "definitions.h"

POINTER init(int xd, int yd, int zd, POINTER V)
{ 
	int i,j,k,mm,a;
	for (i=0; i<xd; i++)
		for (j=0; j<yd; j++)
			for (k=0; k<zd; k++)
				if (Fs(V.s,i,j,k,xd,yd,zd)>-1)
				{
//					Fs(V.rho,i,j,k,xd,yd,zd)=0.;
					for (mm=0;mm<19;mm++)
					{
						Fb(V.f,i,j,k,mm,xd,yd,zd,19)=w[mm];
						Fb(V.ftemp,i,j,k,mm,xd,yd,zd,19)=w[mm];
					}
//					Fb(V.jj,i,j,k,0,xd,yd,zd,3)=Fb(V.jj,i,j,k,1,xd,yd,zd,3)=0.0;
//					Fb(V.jj,i,j,k,2,xd,yd,zd,3)=0.;//0.0577350;
				}
return V;
}