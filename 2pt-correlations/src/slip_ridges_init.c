#include "definitions.h"

POINTER slip_ridges_init(int xd, int yd, int zd, POINTER V, PDATA pd)
{
	int i,j,k,a,b,extx;
	int gap, width;
	int slip=-2,solid=-1;
	int st;
	MPI_Status status;

	gap=gplus;
	width=wplus;
	st=gst;
	extx=xd;

	for (k=0;k<zd;k++)
		for (j=0;j<yd;j++)
			for (i=0;i<xd;i++)
				Fs(V.s,i,j,k,xd,yd,zd)=solid;
	for (k=0;k<zd;k++)
	{
		for (j=st;j<yd;j+=(gap+width))
		{
			for (a=0;a<gap;a++)
				for (i=0;i<extx;i++)
					Fs(V.s,i,(j+a),k,xd,yd,zd)=slip;
		}
	}
	return V;
}