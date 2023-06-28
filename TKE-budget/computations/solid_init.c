#include "definitions.h"

POINTER solid_init(int xd, int yd, int zd, POINTER V)
{
	int i,j,k;
	int solid=-1;
	
	for (i=0;i<xd;i++)
		for (j=0;j<yd;j++)
			for (k=0;k<zd;k++)
				Fs(V.s,i,j,k,xd,yd,zd)=solid;
	return V;
}