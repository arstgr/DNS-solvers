#include "definitions.h"

POINTER reading(int xd, int yd, int zd, POINTER V, PDATA pd, long ts, DP Fx, DP ubulk)
{
	int i,j,k,nn,mm;
	int a,q,b,kk;
	FILE *sv;
	int exty,extx;
	char fn[40];
	long cntr;
	double ux,uy,uz,dens;
	DP rl,rr,dq,kb;
	const int L[5]={2,8,9,13,14};
	const int R[5]={1,7,10,11,12};
	double *af,*sf;
	MPI_File fh;
    MPI_Status status;
    MPI_Offset offset;

	
	cntr=ts;
	if (pd.myrank==0)
		fprintf(stderr,"TRYING TO READ DATA FROM FILES   Fx=%.14f\n",Fx);
	sprintf(fn,"../vel.%.4d",(int)cntr);
        fprintf(stderr,"reading data from %s\n",fn);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

//	offset=b*(xd-1)*3*yd*zd*sizeof(double);
	offset=b*xd*3*yd*zd*sizeof(double);
    MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.vel+2*(yd*zd*3)),((xd-1)*yd*zd*3),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);	
//	offset=(pd.left*(xd-1)*3*yd*zd+(xd-3)*3*yd*zd)*sizeof(double);
	offset=(pd.left*xd*3*yd*zd+(xd-3)*3*yd*zd)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.vel),(yd*zd*3),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);	
//	offset=(pd.left*(xd-1)*3*yd*zd+(xd-2)*3*yd*zd)*sizeof(double);
	offset=(pd.left*xd*3*yd*zd+(xd-2)*3*yd*zd)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.vel+yd*zd*3),(yd*zd*3),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);
//	offset=(pd.right*(xd-1)*3*yd*zd)*sizeof(double);
	offset=(pd.right*xd*3*yd*zd)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.vel+(xd+1)*yd*zd*3),(yd*zd*3),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);
//	offset=(pd.right*(xd-1)*3*yd*zd+yd*zd*3)*sizeof(double);
	offset=(pd.right*xd*3*yd*zd+3*yd*zd)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.vel+(xd+2)*yd*zd*3),(yd*zd*3),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);

	sprintf(fn,"../den.%.4d",(int)cntr);
        fprintf(stderr,"reading data from %s\n",fn);
	MPI_File_open(MPI_COMM_WORLD,fn,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	b=pd.myrank;

//	offset=b*(xd-1)*yd*zd*sizeof(double);
	offset=b*xd*yd*zd*sizeof(double);
    MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.den+2*(yd*zd)),((xd-1)*yd*zd),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);	
//	offset=(pd.left*(xd-1)*yd*zd+(xd-3)*yd*zd)*sizeof(double);
	offset=(pd.left*xd*yd*zd+(xd-3)*yd*zd)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.den),(yd*zd),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);	
//	offset=(pd.left*(xd-1)*yd*zd+(xd-2)*yd*zd)*sizeof(double);
	offset=(pd.left*xd*yd*zd+(xd-2)*yd*zd)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.den+yd*zd),(yd*zd),MPI_DOUBLE,&status);

	MPI_Barrier(MPI_COMM_WORLD);
//	offset=(pd.right*(xd-1)*yd*zd)*sizeof(double);
	offset=(pd.right*xd*yd*zd)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.den+(xd+1)*yd*zd),(yd*zd),MPI_DOUBLE,&status);
	
	MPI_Barrier(MPI_COMM_WORLD);
//	offset=(pd.right*(xd-1)*yd*zd+yd*zd)*sizeof(double);
	offset=(pd.right*xd*yd*zd+yd*zd)*sizeof(double);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_read(fh,(V.den+(xd+2)*yd*zd),(yd*zd),MPI_DOUBLE,&status);
	
	MPI_File_close(&fh);

return V;
}