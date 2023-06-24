#include "definitions.h"

POINTER mean_subtract(int xd,int yd,int zd,POINTER V,PDATA pd,long tst,DP Fx)
{
	DP nu,tau;
	int i,j,k,a,q,b,c;
	FILE *sv;
	int jp,jn,in,ip;
	int jpp,jnn,inn,ipp;
	int num=0,temp=0,upn=0;
	int extx,exty;
	char fn[50],fh[100],fht[100];
	DP uAVG[yd][zd],vAVG[yd][zd],wAVG[yd][zd];
	DP uAVGm[gplus+wplus][zd],vAVGm[gplus+wplus][zd],wAVGm[gplus+wplus][zd];
	DP uAVG2[yd][(zd+1)/2],vAVG2[yd][(zd+1)/2],wAVG2[yd][(zd+1)/2];
	DP *UAVG,*VAVG,*WAVG,*DUMMY,*DUMMYYZ;
	MPI_Status status;
	DP UA[zd],VA[zd],WA[zd];
		

	extx=xd;
	exty=yd;
	upn=1+(zd-1)/2;
	num=extx*exty*pd.numproc;

	tau=(3./Rebs)*ub*(zd-2.)+0.5;
	nu=(1./6.)*(2.*tau-1.);  
	
	fprintf(stderr,"starting mean subtract i'm %d\n",pd.myrank);
	
	DUMMYYZ=(DP *)calloc(yd*zd,sizeof(DP));
	
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			uAVG[j][k]=0.;
			vAVG[j][k]=0.;
			wAVG[j][k]=0.;
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

	for (a=0;a<(gplus+wplus);a++)
	{
		for (k=0;k<zd;k++)
		{
			uAVGm[a][k]=0.;
			vAVGm[a][k]=0.;
			wAVGm[a][k]=0.;
		}
	}

	for (j=0;j<yd;j++)
	{
		for (k=0;k<upn;k++)
		{
			uAVG2[j][k]=0.;
                        vAVG2[j][k]=0.;
                        wAVG2[j][k]=0.;
		}
	}
	
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
		}
	}

	
	for (j=0;j<yd;j++)
	{
		for (k=1;k<upn;k++)
		{
			uAVG2[j][k]=0.5*((uAVG[j][k])+(uAVG[j][zd-1-k]));
			vAVG2[j][k]=0.5*((vAVG[j][k])+(vAVG[j][zd-1-k]));
			wAVG2[j][k]=0.5*((wAVG[j][k])+(-wAVG[j][zd-1-k]));
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
		}
	}
/*
	for (i=0;i<extx;i++)
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
*/
	for (k=0;k<zd;k++)
	{
		UA[k]=0.;
		for (j=0;j<exty;j++)
			UA[k]+=uAVG[j][k];
		UA[k]/=((double)exty);

		VA[k]=0.;
                for (j=0;j<exty;j++)
                        VA[k]+=vAVG[j][k];
                VA[k]/=((double)exty);

		WA[k]=0.;
                for (j=0;j<exty;j++)
                        WA[k]+=wAVG[j][k];
                WA[k]/=((double)exty);
	}


	for (i=0;i<extx;i++)
        {
                for (j=0;j<exty;j++)
                {
                        for (k=1;k<(zd-1);k++)
                        {
                                Fb(V.vel,i,j,k,0,xd,yd,zd,3) -= UA[k];
                                Fb(V.vel,i,j,k,1,xd,yd,zd,3) -= VA[k];
                                Fb(V.vel,i,j,k,2,xd,yd,zd,3) -= WA[k];
                        }
                }
        }

	if (pd.myrank==0)
	{
		sprintf(fn,"uavg.%.3d",((int)tst));
		sv=fopen(fn,"wb");
		fwrite(UA,sizeof(double),zd,sv);
		fclose(sv);

/*		sprintf(fn,"vavg.%.3d",((int)tst));
                sv=fopen(fn,"wb");
                fwrite(VA,sizeof(double),zd,sv);
                fclose(sv);

		sprintf(fn,"wavg.%.3d",((int)tst));
                sv=fopen(fn,"wb");
                fwrite(WA,sizeof(double),zd,sv);
                fclose(sv);       */
	}
	
	free(DUMMYYZ);
	fprintf(stderr,"end of mean subtract i'm %d\n",pd.myrank);
  
  return V;
}