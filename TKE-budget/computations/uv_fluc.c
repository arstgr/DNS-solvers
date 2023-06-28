#include "definitions.h"

POINTER uv_fluc(int xd, int yd, int zd, POINTER V, PDATA pd, long ts, DP Fx, int position)
{
	int i,j,k,a,q,c,upn;
	DP uAVG[yd][zd],vAVG[yd][zd],wAVG[yd][zd],dAVG[yd][zd];
	DP uAVGm[gplus+wplus][zd],vAVGm[gplus+wplus][zd],wAVGm[gplus+wplus][zd],dAVGm[gplus+wplus][zd];
	DP uAVG2[yd][(zd+1)/2],vAVG2[yd][(zd+1)/2],wAVG2[yd][(zd+1)/2],dAVG2[yd][(zd+1)/2];
	double *DUMMYYZ;
	double ufl[xd-1][yd],vfl[xd-1][yd];
	char fl[30];
	MPI_File fh;
	MPI_Offset offset;
	MPI_Status status;
	int extx=xd-1,exty=yd;

	upn=1+(zd-1)/2;
	offset=pd.myrank*(xd-1)*yd*sizeof(double);
	DUMMYYZ=(DP *)calloc(yd*zd,sizeof(DP));
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			uAVG[j][k]=0.;
			vAVG[j][k]=0.;
			wAVG[j][k]=0.;
			dAVG[j][k]=0.;
		}
	}
		
		
	for (j=0;j<exty;j++)
	{
		for (k=1;k<(zd-1);k++)
		{
			for (i=0;i<extx;i++)
			{
				uAVG[j][k]+=Fb(V.vel,i,j,k,0,xd,yd,zd,3);
				vAVG[j][k]+=Fb(V.vel,i,j,k,1,xd,yd,zd,3);
				wAVG[j][k]+=Fb(V.vel,i,j,k,2,xd,yd,zd,3);
				dAVG[j][k]+=(1./3.)*Fb(V.den,i,j,k,0,xd,yd,zd,1);
			}
		}
	}
		
	for (j=0;j<exty;j++)
	{
		for (k=1;k<(zd-1);k++)
		{
			uAVG[j][k]/=(extx*pd.numproc);
			vAVG[j][k]/=(extx*pd.numproc);
			wAVG[j][k]/=(extx*pd.numproc);
			dAVG[j][k]/=(extx*pd.numproc);
		}
	}
	
	MPI_Reduce(&uAVG[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uAVG,DUMMYYZ,yd*zd*sizeof(DP));	
	
	MPI_Reduce(&vAVG[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vAVG,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&wAVG[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wAVG,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&dAVG[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(dAVG,DUMMYYZ,yd*zd*sizeof(DP));

	for (a=0;a<(gplus+wplus);a++)
	{
		for (k=0;k<zd;k++)
		{
			uAVGm[a][k]=0.;
			vAVGm[a][k]=0.;
			wAVGm[a][k]=0.;
			dAVGm[a][k]=0.;
		}
	}

	for (j=0;j<yd;j++)
	{
		for (k=0;k<upn;k++)
		{
			uAVG2[j][k]=0.;
                        vAVG2[j][k]=0.;
                        wAVG2[j][k]=0.;
                        dAVG2[j][k]=0.;
		}
	}

// ****************************************************************************************

	for (j=gst;j<yd;j+=(gplus+wplus))
	{
		for (a=0;(a<(gplus+wplus))&&((a+j)<yd);a++)
		{
			for (k=0;k<zd;k++)
			{
				c=(a+j-gst)%(gplus+wplus);
				uAVGm[c][k]+=uAVG[a+j][k];
				vAVGm[c][k]+=vAVG[a+j][k];
				wAVGm[c][k]+=wAVG[a+j][k];
				dAVGm[c][k]+=dAVG[a+j][k];
			}
		}
	}

	for (a=0;a<gst;a++)
	{
		for (k=0;k<zd;k++)
		{
			c=a+gst+gplus;
			uAVGm[c][k]+=uAVG[a][k];
			vAVGm[c][k]+=vAVG[a][k];
			wAVGm[c][k]+=wAVG[a][k];
			dAVGm[c][k]+=dAVG[a][k];
		}
	}

	c=yd/(gplus+wplus);
	for (a=0;a<(gplus+wplus);a++)
	{
		for (k=0;k<zd;k++)
		{
			uAVGm[a][k]/=(double)c;
			vAVGm[a][k]/=(double)c;
			wAVGm[a][k]/=(double)c;
			dAVGm[a][k]/=(double)c;
		}
	}

	for (j=gst;j<yd;j+=(gplus+wplus))
	{
		for (a=0;(a<(gplus+wplus))&&((a+j)<yd);a++)
		{
			for (k=0;k<zd;k++)
			{
				c=(a+j-gst)%(gplus+wplus);
				uAVG[j+a][k]=uAVGm[c][k];
				vAVG[j+a][k]=vAVGm[c][k];
				wAVG[j+a][k]=wAVGm[c][k];
				dAVG[j+a][k]=dAVGm[c][k];
			}
		}
	}

	for (a=0;a<gst;a++)
	{
		for (k=0;k<zd;k++)
		{
			c=a+gst+gplus;
			uAVG[a][k]=uAVGm[c][k];
			vAVG[a][k]=vAVGm[c][k];
			wAVG[a][k]=wAVGm[c][k];
			dAVG[a][k]=dAVGm[c][k];
		}
	}

	for (j=0;j<yd;j++)
	{
		for (k=1;k<upn;k++)
		{
			uAVG2[j][k]=0.5*((uAVG[j][k])+(uAVG[j][zd-1-k]));
			vAVG2[j][k]=0.5*((vAVG[j][k])+(vAVG[j][zd-1-k]));
			wAVG2[j][k]=0.5*((wAVG[j][k])+(-wAVG[j][zd-1-k]));
			dAVG2[j][k]=0.5*((dAVG[j][k])+(dAVG[j][zd-1-k]));
		}
	}

	for (j=0;j<yd;j++)
	{
		for (k=1;k<upn;k++)
		{
			uAVG[j][k]=uAVG2[j][k];
			uAVG[j][zd-1-k]=uAVG2[j][k];
			vAVG[j][k]=vAVG2[j][k];
       			vAVG[j][zd-1-k]=vAVG2[j][k];
			wAVG[j][k]=wAVG2[j][k];
			wAVG[j][zd-1-k]=-wAVG2[j][k];
			dAVG[j][k]=dAVG2[j][k];
			dAVG[j][zd-1-k]=dAVG2[j][k];
		}
	}

//**************************************************************************************************************

	for (i=-2;i<(extx+2);i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=1;k<(zd-1);k++)
			{
				Fb(V.vel,i,j,k,0,xd,yd,zd,3) -= uAVG[j][k];
				Fb(V.vel,i,j,k,1,xd,yd,zd,3) -= vAVG[j][k];
				Fb(V.vel,i,j,k,2,xd,yd,zd,3) -= wAVG[j][k];
			}
		}
	}
	
	k=position;
	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			ufl[i][j] = Fb(V.vel,i,j,k,0,xd,yd,zd,3);
			vfl[i][j] = Fb(V.vel,i,j,k,1,xd,yd,zd,3);
		}
	}

	sprintf(fl,"ufl4.%.4d",((int)ts));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_write_all(fh,&ufl[0][0],(xd-1)*yd,MPI_DOUBLE,&status);
	MPI_File_close(&fh);

	sprintf(fl,"vfl4.%.4d",((int)ts));
	MPI_File_open(MPI_COMM_WORLD,fl,MPI_MODE_CREATE | MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_set_view(fh,offset,MPI_DOUBLE,MPI_DOUBLE,"native",MPI_INFO_NULL);
	MPI_File_write_all(fh,&vfl[0][0],(xd-1)*yd,MPI_DOUBLE,&status);
	MPI_File_close(&fh);

	free(DUMMYYZ);
return V;
}