#include "definitions.h"

POINTER solid_init(int xd, int yd, int zd, POINTER V)
{
	int i,j,k;
	int solid=-1,periodic=1,fluid=0;
	
	for (i=0;i<xd;i++)
		for (j=0;j<yd;j++)
			for (k=0;k<zd;k++)
				Fs(V.s,i,j,k,xd,yd,zd)=fluid;
	
	k=0;
	for (i=0;i<xd;i++)
		for (j=0;j<yd;j++)
			Fs(V.s,i,j,k,xd,yd,zd)=solid;
	k=zd-1;
	for (i=0;i<xd;i++)
		for (j=0;j<yd;j++)
			Fs(V.s,i,j,k,xd,yd,zd)=solid;
	
	i=0;
	for (j=0;j<yd;j++)
		for (k=1;k<(zd-1);k++)
			Fs(V.s,i,j,k,xd,yd,zd)=periodic;
	i=xd-1;
	for (j=0;j<yd;j++)
		for (k=1;k<(zd-1);k++)
			Fs(V.s,i,j,k,xd,yd,zd)=periodic;
	
	j=0;
	for (i=0;i<xd;i++)
		for (k=1;k<(zd-1);k++)
			Fs(V.s,i,j,k,xd,yd,zd)=periodic;
	j=yd-1;
	for (i=0;i<xd;i++)
		for (k=1;k<(zd-1);k++)
			Fs(V.s,i,j,k,xd,yd,zd)=periodic;
	return V;
}