#include "definitions.h"

POINTER reading(int xd, int yd, int zd, POINTER V, PDATA pd, long ts, DP Fx, DP ubulk)
{
	int i,j,k,nn,mm;
	int a,q,b,kk;
	FILE *sv;
	int exty,extx;
	char fn[40];
	long cntr;
	MPI_File fh;
    MPI_Status status;
    MPI_Offset offset;

	cntr=ts;
	if (pd.myrank==0)
		fprintf(stderr,"TRYING TO READ DATA FROM FILES, I'm root\n");
	sprintf(fn,"../vel.%.4d",cntr);
        fprintf(stderr,"reading data from %s myrank is%d\n",fn,pd.myrank);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

	offset=b*xd*3*yd*zd*sizeof(double);
//	offset=b*(xd+1)*3*yd*zd*sizeof(double);
    MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,V.vel,(xd*yd*zd*3),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);

	fprintf(stderr,"reading finished i'm %d\n",pd.myrank);
return V;
}