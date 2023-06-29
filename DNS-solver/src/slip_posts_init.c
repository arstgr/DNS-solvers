#include "definitions.h"

POINTER slip_posts_init(int xd, int yd, int zd, POINTER V, PDATA pd)
{
        int i,j,k,a,b,extx,exty,extn,t1,t2;
        int gx,gy, wx,wy;
        int slip=-2;
        MPI_Status status;
        int solid=-1,periodic=1,fluid=0;

        extx=xd-1;
        extn=(xd-1)*pd.numproc;

        gx=2;
        wx=2;

        gy=2;
        wy=2;
	
	exty =yd-1;

        MPI_Barrier(MPI_COMM_WORLD);
        k=0;
        for (j=1;j<exty;j+=(gy+wy))
        {
                for (a=0;a<gy;a++)
                {
                        for (i=1;i<extn;i+=(gx+wx))
                        {
                                for (b=0;b<gx;b++)
                                {
                                        t1=i-pd.myrank*extx;
                                        if ((t1>=0)&&(t1<xd))
                                        {
                                                if (((t1+b)<xd)&&((j+a)<yd))
                                                        Fs(V.s,(t1+b),(j+a),k,xd,yd,zd)=slip;
                                        }
                                }
                        }
                }
        }
	j=yd-1;
	for (i=0;i<extx;i++)
		Fs(V.s,i,j,k,xd,yd,zd)=Fs(V.s,i,(j-1),k,xd,yd,zd);
	
        MPI_Barrier(MPI_COMM_WORLD);
        k=zd-1;
        for (j=1;j<exty;j+=(gy+wy))
        {
                for (a=0;a<gy;a++)
                {
                        for (i=1;i<extn;i+=(gx+wx))
                        {
                                for (b=0;b<gx;b++)
                                {
                                        t1=i-pd.myrank*extx;
                                        if ((t1>=0)&&(t1<xd))
                                        {
                                                if (((t1+b)<xd)&&((j+a)<yd))
                                                        Fs(V.s,(t1+b),(j+a),k,xd,yd,zd)=slip;
                                        }
                                }
                        }
                }
        }
        	
	j=yd-1;
	for (i=0;i<extx;i++)
		Fs(V.s,i,j,k,xd,yd,zd)=Fs(V.s,i,(j-1),k,xd,yd,zd);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Sendrecv(V.s,(yd*zd),MPI_INT,pd.left,1169,(V.s+(xd-1)*yd*zd),(yd*zd),MPI_INT,pd.right,1169,MPI_COMM_WORLD,&status);

        return V;
}