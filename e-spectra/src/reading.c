#include "definitions.h"

POINTER reading(int xd, int yd, int zd, POINTER V, PDATA pd, int ts)
{
	char fn[40];
	int cntr,b;
	MPI_File fh;
	MPI_Status status;
	MPI_Offset offset,off1,off2;

	cntr=ts;
	if (pd.myrank==0)
		fprintf(stderr,"TRYING TO READ DATA FROM FILES\n");
	sprintf(fn,"../vel.%.4d",cntr);
	fprintf(stderr,"reading data from %s\n",fn);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

//	offset=b*(xd-1)*3*yd*zd*sizeof(double);
	off1 = (MPI_Offset)(zd);
	off2 = (MPI_Offset)(yd)*off1;
	offset=((MPI_Offset)b)*((MPI_Offset)(xd-1))*((MPI_Offset)3)*off2*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.vel),((xd-1)*yd*zd*3),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);

  return V;
}