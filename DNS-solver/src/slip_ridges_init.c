#include "definitions.h"

POINTER slip_ridges_init(int xd, int yd, int zd, POINTER V, PDATA pd, MPI_Comm cart_grid)
{
	int i,j,k,a,b,extx,exty;
	int gap, width,halfw;
	int slip=-2;
	MPI_Status status;
	int lim1,lim2,coord[2],jstar;

	gap=2;
	width=2;
	halfw=width/2;
	extx=xd;
	exty = (yd-1)*YDIM;
	
	MPI_Cart_coords(cart_grid, pd.myrank, 2, coord);
	lim1 = coord[1]*(yd-1);
	lim2 = (coord[1]+1)*(yd-1);

	k=0;
	for (j=halfw;j<exty;j+=(gap+width))
	{
		if ((j>=lim1)&&(j<lim2))
		{
			jstar = j%(yd-1);
			for (a=0;a<gap;a++)
				for (i=0;i<extx;i++)
					Fs(V.s,i,(jstar+a),k,xd,yd,zd)=slip;
		}
	}
	k=0;
	j=yd-1;
	for (i=0;i<extx;i++)
		Fs(V.s,i,j,k,xd,yd,zd)=Fs(V.s,i,(j-1),k,xd,yd,zd);
	
	k=zd-1;
	for (i=0;i<extx;i++)
	  for (j=0;j<yd;j++)
	    Fs(V.s,i,j,k,xd,yd,zd)=Fs(V.s,i,j,0,xd,yd,zd);

	return V;
}
