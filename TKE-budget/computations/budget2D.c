#include "definitions.h"

POINTER budget2D(int xd, int yd, int zd, POINTER V, PDATA pd, long ts, DP Fx)
{
	int i,j,k,a,q,c;
	DP ufl,vfl,wfl,dfl,dens;
	DP uAVG[yd][zd],vAVG[yd][zd],wAVG[yd][zd],dAVG[yd][zd];
	DP uAVGm[gplus+wplus][zd],vAVGm[gplus+wplus][zd],wAVGm[gplus+wplus][zd],dAVGm[gplus+wplus][zd];
	DP uAVG2[yd][(zd+1)/2],vAVG2[yd][(zd+1)/2],wAVG2[yd][(zd+1)/2],dAVG2[yd][(zd+1)/2];
	DP *UAVG,*VAVG,*WAVG,*DUMMY,*DUMMYYZ;
	DP *UTIN,*VTIN,*WTIN,*DTIN;
	double dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz;
	double *uu,*uv,*uw,*vv,*vw,*ww,*dd;
	DP *DUMMYXZ,*DAVGXZ,*DAVG;
	DP siu,siv,siw,sip;
	double *SKEWu,*SKEWv,*SKEWw,*SKEWp,*SKEWuw,*SKEWuv,*SKEWvw;
	double *KURTu,*KURTv,*KURTw,*KURTp,*KURTuw,*KURTuv,*KURTvw;
	double ox,oy,oz;
	double u,v,w;
	double *wx,*wy,*wz;
	double *UW;
	double *KE,KEt;
	double *Pii,*Puu,*Pvv,*Pww,*Puw,*Puv,*Pvw;
	double *Tiir,*Tiipai,*Tiis,*Epii;
	double *Paiuu,*Paivv,*Paiww,*Paiuw,*Paiuv,*Paivw,*Epuu,*Epvv,*Epww,*Epuw,*Epuv,*Epvw;
	double *Tiirt,*Tuur,*Tvvr,*Twwr,*Tuwr,*Tuvr,*Tvwr;
	double *Tuupai,*Tvvpai,*Twwpai,*Tuwpai,*Tuvpai,*Tvwpai;
	double *Tuus,*Tvvs,*Twws,*Tuws,*Tuvs,*Tvws;
	double tmp;
	double Pii2d[yd][zd],Tiipai2d[yd][zd],Tiis2d[yd][zd],Tiir2d[yd][zd],Epii2d[yd][zd];
	double Puu2d[yd][zd],Tuupai2d[yd][zd],Tuus2d[yd][zd],Tuur2d[yd][zd],Epuu2d[yd][zd];
	double Pvv2d[yd][zd],Tvvpai2d[yd][zd],Tvvs2d[yd][zd],Tvvr2d[yd][zd],Epvv2d[yd][zd];
	double Pww2d[yd][zd],Twwpai2d[yd][zd],Twws2d[yd][zd],Twwr2d[yd][zd],Epww2d[yd][zd];
	double Puv2d[yd][zd],Tuvpai2d[yd][zd],Tuvs2d[yd][zd],Tuvr2d[yd][zd],Epuv2d[yd][zd];
	double Puw2d[yd][zd],Tuwpai2d[yd][zd],Tuws2d[yd][zd],Tuwr2d[yd][zd],Epuw2d[yd][zd];
	double Pvw2d[yd][zd],Tvwpai2d[yd][zd],Tvws2d[yd][zd],Tvwr2d[yd][zd],Epvw2d[yd][zd];
	double KE2d[yd][zd],dUdy2d[yd][zd],dUdz2d[yd][zd],dVdy2d[yd][zd],dVdz2d[yd][zd],dWdy2d[yd][zd],dWdz2d[yd][zd];
	double uu2d[yd][zd],uv2d[yd][zd],uw2d[yd][zd],vv2d[yd][zd],vw2d[yd][zd],ww2d[yd][zd];
	double pu2d[yd][zd],pv2d[yd][zd],pw2d[yd][zd];
	double Paiuu2d[yd][zd],Paivv2d[yd][zd],Paiww2d[yd][zd],Paiuv2d[yd][zd],Paiuw2d[yd][zd],Paivw2d[yd][zd];
	double pdudx2d[yd][zd],pdvdy2d[yd][zd],pdwdz2d[yd][zd],pdwdy2d[yd][zd];
	double pdvdx2d[yd][zd],pdudy2d[yd][zd],pdudz2d[yd][zd],pdwdx2d[yd][zd],pdvdz2d[yd][zd];
	double usuv2d[yd][zd],vsvv2d[yd][zd],wswv2d[yd][zd],wsuv2d[yd][zd],usvw2d[yd][zd];
	double usuw2d[yd][zd],vsvw2d[yd][zd],wsww2d[yd][zd],vsuw2d[yd][zd],wsvv2d[yd][zd];
	double usvv2d[yd][zd],vsuv2d[yd][zd],uswv2d[yd][zd],usww2d[yd][zd],wsuw2d[yd][zd];
	double vsww2d[yd][zd];
	double vk2d[yd][zd],wk2d[yd][zd];
	double VdKdy2d[yd][zd],WdKdz2d[yd][zd];
	double vke2d[yd][zd],wke2d[yd][zd];
	double uuv2d[yd][zd],uuw2d[yd][zd],vvv2d[yd][zd],vvw2d[yd][zd],wwv2d[yd][zd],www2d[yd][zd];
	double uvv2d[yd][zd],uvw2d[yd][zd],uwv2d[yd][zd],uww2d[yd][zd],vwv2d[yd][zd],vww2d[yd][zd];
	double *VdKdy,*WdKdz;
	double *Vduudy,*Wduudz,*Vdvvdy,*Wdvvdz,*Vdwwdy,*Wdwwdz,*Vduvdy,*Wduvdz,*Vduwdy,*Wduwdz,*Vdvwdy,*Wdvwdz;
	double Vduudy2d[yd][zd],Wduudz2d[yd][zd],Vdvvdy2d[yd][zd],Wdvvdz2d[yd][zd],Vdwwdy2d[yd][zd],Wdwwdz2d[yd][zd];
	double Vduvdy2d[yd][zd],Vduwdy2d[yd][zd],Vdvwdy2d[yd][zd],Wduvdz2d[yd][zd],Wduwdz2d[yd][zd],Wdvwdz2d[yd][zd];
	double P1[yd][zd],P2[yd][zd],P3[yd][zd],P4[yd][zd],P5[yd][zd],P6[yd][zd],P7[yd][zd],P8[yd][zd],P9[yd][zd];
	double P10[yd][zd],P11[yd][zd],P12[yd][zd],P13[yd][zd],P14[yd][zd],P15[yd][zd],P16[yd][zd],P17[yd][zd],P18[yd][zd];
	double p1[zd],p2[zd],p3[zd],p4[zd],p5[zd],p6[zd],p7[zd],p8[zd],p9[zd],p10[zd],p11[zd],p12[zd],p13[zd],p14[zd],p15[zd];
	double p16[zd],p17[zd],p18[zd];
	double UW2d[yd][zd];
	double OmegaX[yd][zd];
	double omegax2d[yd][zd],omegay2d[yd][zd],omegaz2d[yd][zd];
	double dudzusuw[yd][zd];
	double dd2d[yd][zd];
	DP nu,tau;
	FILE *sv;
	int jp,jn,in,ip;
	int jpp,jnn,inn,ipp;
	int num=0,temp=0,upn;
	int extx,exty;
	char fn[50],fh[100],fht[100];
	MPI_Status status;
	
	double quadrant[4][zd],dquadrant[4][zd],quadrant2d[4][yd][zd];
	int nquadrant[4][zd];
	double *DUMMYq,*DUMMYq2d;
	int *iDUMMYq;

	extx=xd-1;
	exty=yd;
	upn=1+(zd-1)/2;
	num=extx*exty*pd.numproc;

	tau=(3./Rebs)*ub*(zd-2.)+0.5;
	nu=(1./6.)*(2.*tau-1.);
	
	UAVG=(DP *)calloc(zd,sizeof(DP));
	VAVG=(DP *)calloc(zd,sizeof(DP));
	WAVG=(DP *)calloc(zd,sizeof(DP));
	DAVG=(DP *)calloc(zd,sizeof(DP));
	DUMMY=(DP *)calloc(zd,sizeof(DP));
	UW=(DP *)calloc(zd,sizeof(DP));	

	DUMMYYZ=(DP *)calloc(yd*zd,sizeof(DP));
	DUMMYq=(DP *)calloc(4*zd,sizeof(DP));
	iDUMMYq=(int *)calloc(4*zd,sizeof(DP));
	DUMMYq2d=(DP *)calloc(4*yd*zd,sizeof(DP));

	UTIN=(DP *)calloc(zd,sizeof(DP));
	VTIN=(DP *)calloc(zd,sizeof(DP));
	WTIN=(DP *)calloc(zd,sizeof(DP));
	DTIN=(DP *)calloc(zd,sizeof(DP));
	
	uu=(DP *)calloc(zd,sizeof(DP));
	uv=(DP *)calloc(zd,sizeof(DP));
	uw=(DP *)calloc(zd,sizeof(DP));
	vv=(DP *)calloc(zd,sizeof(DP));
	vw=(DP *)calloc(zd,sizeof(DP));
	ww=(DP *)calloc(zd,sizeof(DP));
	dd=(DP *)calloc(zd,sizeof(DP));
	
	wx=(DP *)calloc(zd,sizeof(DP));
	wy=(DP *)calloc(zd,sizeof(DP));
	wz=(DP *)calloc(zd,sizeof(DP));
	
	SKEWu=(DP *)calloc(zd,sizeof(DP));
	SKEWv=(DP *)calloc(zd,sizeof(DP));
	SKEWw=(DP *)calloc(zd,sizeof(DP));
	SKEWp=(DP *)calloc(zd,sizeof(DP));
	SKEWuw=(DP *)calloc(zd,sizeof(DP));
	SKEWuv=(DP *)calloc(zd,sizeof(DP));
	SKEWvw=(DP *)calloc(zd,sizeof(DP));
	
	KURTu=(DP *)calloc(zd,sizeof(DP));
	KURTv=(DP *)calloc(zd,sizeof(DP));
	KURTw=(DP *)calloc(zd,sizeof(DP));
	KURTp=(DP *)calloc(zd,sizeof(DP));
	KURTuw=(DP *)calloc(zd,sizeof(DP));
	KURTvw=(DP *)calloc(zd,sizeof(DP));
	KURTuv=(DP *)calloc(zd,sizeof(DP));

	KE=(DP *)calloc(zd,sizeof(DP));

	Paiuu=(DP *)calloc(zd,sizeof(DP));
	Paivv=(DP *)calloc(zd,sizeof(DP));
	Paiww=(DP *)calloc(zd,sizeof(DP));
	Paiuw=(DP *)calloc(zd,sizeof(DP));
	Paiuv=(DP *)calloc(zd,sizeof(DP));
	Paivw=(DP *)calloc(zd,sizeof(DP));

	Tiirt=(DP *)calloc(zd,sizeof(DP));
	Tuur=(DP *)calloc(zd,sizeof(DP));
	Tvvr=(DP *)calloc(zd,sizeof(DP));
	Twwr=(DP *)calloc(zd,sizeof(DP));
	Tuwr=(DP *)calloc(zd,sizeof(DP));
	Tuvr=(DP *)calloc(zd,sizeof(DP));
	Tvwr=(DP *)calloc(zd,sizeof(DP));

	Tuupai=(DP *)calloc(zd,sizeof(DP));
	Tvvpai=(DP *)calloc(zd,sizeof(DP));
	Twwpai=(DP *)calloc(zd,sizeof(DP));
	Tuwpai=(DP *)calloc(zd,sizeof(DP));
	Tuvpai=(DP *)calloc(zd,sizeof(DP));
	Tvwpai=(DP *)calloc(zd,sizeof(DP));

	Tuus=(DP *)calloc(zd,sizeof(DP));
	Tvvs=(DP *)calloc(zd,sizeof(DP));
	Twws=(DP *)calloc(zd,sizeof(DP));
	Tuws=(DP *)calloc(zd,sizeof(DP));
	Tuvs=(DP *)calloc(zd,sizeof(DP));
	Tvws=(DP *)calloc(zd,sizeof(DP));

	Epuu=(DP *)calloc(zd,sizeof(DP));
	Epvv=(DP *)calloc(zd,sizeof(DP));
	Epww=(DP *)calloc(zd,sizeof(DP));
	Epuw=(DP *)calloc(zd,sizeof(DP));
	Epuv=(DP *)calloc(zd,sizeof(DP));
	Epvw=(DP *)calloc(zd,sizeof(DP));

	Pii=(DP *)calloc(zd,sizeof(DP));
	Puu=(DP *)calloc(zd,sizeof(DP));
	Pvv=(DP *)calloc(zd,sizeof(DP));
	Pww=(DP *)calloc(zd,sizeof(DP));
	Puw=(DP *)calloc(zd,sizeof(DP));
	Puv=(DP *)calloc(zd,sizeof(DP));
	Pvw=(DP *)calloc(zd,sizeof(DP));
	
	Tiir=(DP *)calloc(zd,sizeof(DP));
	Tiipai=(DP *)calloc(zd,sizeof(DP));
	Tiis=(DP *)calloc(zd,sizeof(DP));
	Epii=(DP *)calloc(zd,sizeof(DP));
	
	VdKdy=(DP *)calloc(zd,sizeof(DP));
	WdKdz=(DP *)calloc(zd,sizeof(DP));
	
	Vduudy=(DP *)calloc(zd,sizeof(DP));
	Wduudz=(DP *)calloc(zd,sizeof(DP));
	Vdvvdy=(DP *)calloc(zd,sizeof(DP));
	Wdvvdz=(DP *)calloc(zd,sizeof(DP));
	Vdwwdy=(DP *)calloc(zd,sizeof(DP));
	Wdwwdz=(DP *)calloc(zd,sizeof(DP));
	Vduvdy=(DP *)calloc(zd,sizeof(DP));
	Wduvdz=(DP *)calloc(zd,sizeof(DP));
	Vduwdy=(DP *)calloc(zd,sizeof(DP));
	Wduwdz=(DP *)calloc(zd,sizeof(DP));
	Vdvwdy=(DP *)calloc(zd,sizeof(DP));
	Wdvwdz=(DP *)calloc(zd,sizeof(DP));
	
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

	for (j=0;j<exty;j++)
	{			
		for (k=1;k<(zd-1);k++)
		{
			*(UAVG+k) += uAVG[j][k];
			*(VAVG+k) += vAVG[j][k];
			*(WAVG+k) += wAVG[j][k];
			*(DAVG+k) += dAVG[j][k];
			*(UW+k) += uAVG[j][k]*wAVG[j][k];
			UW2d[j][k] = uAVG[j][k]*wAVG[j][k];
		}
	}

		
	for (k=0;k<zd;k++)
	{
		*(UAVG+k) /= (double)(exty);
		*(VAVG+k) /= (double)(exty);
		*(WAVG+k) /= (double)(exty);
		*(DAVG+k) /= (double)(exty);
		*(UW+k) /= (double)(exty);
	}	
			
// *************************************************************************************
/*			
	for (j=0;j<exty;j++)
	{
		for (k=0;k<zd;k++)
		{
//			uAVG[j][k] = (*(UAVG+k));
//			vAVG[j][k] = (*(VAVG+k));
//			wAVG[j][k] = (*(WAVG+k));
			dAVG[j][k] = (*(DAVG+k));
		}
	}
*/	
// *************************************************************************************
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
	
	for (j=0;j<exty;j++)
	{
		for (k=1;k<(zd-1);k++)
		{
			KE2d[j][k] = 0.;
			uu2d[j][k] = 0.;
			uv2d[j][k] = 0.;
			uw2d[j][k] = 0.;
			vv2d[j][k] = 0.;
			vw2d[j][k] = 0.;
			ww2d[j][k] = 0.;
			
			usuv2d[j][k] = 0.;
			usuw2d[j][k] = 0.;
			vsvv2d[j][k] = 0.;
			wswv2d[j][k] = 0.;
			usvw2d[j][k] = 0.;
			vsvw2d[j][k] = 0.;
			vsww2d[j][k] = 0.;
			wsww2d[j][k] = 0.;
			wsuv2d[j][k] = 0.;
			usww2d[j][k] = 0.;
			vsuv2d[j][k] = 0.;
			wsvv2d[j][k] = 0.;
			wsuw2d[j][k] = 0.;
			usvv2d[j][k] = 0.;
			vsuw2d[j][k] = 0.;
			
			pdudx2d[j][k] = 0.;
			pdvdy2d[j][k] = 0.;
			pdwdz2d[j][k] = 0.;
			pdwdy2d[j][k] = 0.;
			pdvdx2d[j][k] = 0.;
			pdudy2d[j][k] = 0.;
			pdudz2d[j][k] = 0.;
			pdwdx2d[j][k] = 0.;
			pdvdz2d[j][k] = 0.;
				
			vke2d[j][k] = 0.;
			wke2d[j][k] = 0.;
			
			uuv2d[j][k] = 0.;
			uuw2d[j][k] = 0.;
			vvv2d[j][k] = 0.;
			vvw2d[j][k] = 0.;
			wwv2d[j][k] = 0.;
			www2d[j][k] = 0.;
			uvv2d[j][k] = 0.;
			uvw2d[j][k] = 0.;
			uwv2d[j][k] = 0.;
			uww2d[j][k] = 0.;
			vwv2d[j][k] = 0.;
			vww2d[j][k] = 0.;
				
			Epii2d[j][k] = 0.;
			Epuu2d[j][k] = 0.;
			Epvv2d[j][k] = 0.;
			Epww2d[j][k] = 0.;
			Epuw2d[j][k] = 0.;
			Epuv2d[j][k] = 0.;
			Epvw2d[j][k] = 0.;
			
			pu2d[j][k] = 0.;
			pv2d[j][k] = 0.;
			pw2d[j][k] = 0.;
			
			OmegaX[j][k] = 0.;
			
			omegax2d[j][k] = 0.;
			omegay2d[j][k] = 0.;
			omegaz2d[j][k] = 0.;

			dudzusuw[j][k] = 0.;

			dd2d[j][k] = 0.;
			
			quadrant2d[0][j][k] = 0.;
			quadrant2d[1][j][k] = 0.;
			quadrant2d[2][j][k] = 0.;
			quadrant2d[3][j][k] = 0.;
		}
	}
	for (k=0;k<zd;k++)
	{
		quadrant[0][k] = 0.;
		quadrant[1][k] = 0.;
		quadrant[2][k] = 0.;
		quadrant[3][k] = 0.;
		
		nquadrant[0][k] = 0;
		nquadrant[1][k] = 0;
		nquadrant[2][k] = 0;
		nquadrant[3][k] = 0;
	}
	
	for (i=0;i<extx;i++)
	{
		for (j=0;j<exty;j++)
		{
			for (k=1;k<(zd-1);k++)
			{
				
				u=Fb(V.vel,i,j,k,0,xd,yd,zd,3);
				v=Fb(V.vel,i,j,k,1,xd,yd,zd,3);
				w=Fb(V.vel,i,j,k,2,xd,yd,zd,3);

				dens = Fb(V.den,i,j,k,0,xd,yd,zd,1);
				dfl = (1./3.)*dens-dAVG[j][k];

				*(uu+k) += u*u;
				*(uv+k) += u*v;
				*(uw+k) += u*w;
				*(vv+k) += v*v;
				*(vw+k) += v*w;
				*(ww+k) += w*w;
				
				*(SKEWu+k) += u*u*u;
				*(SKEWv+k) += v*v*v;
				*(SKEWw+k) += w*w*w;
				*(SKEWp+k) += dfl*dfl*dfl;
				*(SKEWuw+k) += u*u*u*w*w*w;
				*(SKEWuv+k) += u*u*u*v*v*v;
				*(SKEWvw+k) += v*v*v*w*w*w;
				
				*(KURTu+k) += u*u*u*u;
				*(KURTv+k) += v*v*v*v;
				*(KURTw+k) += w*w*w*w;
				*(KURTp+k) += dfl*dfl*dfl*dfl;
				*(KURTuw+k) += u*u*u*u*w*w*w*w;
				*(KURTvw+k) += v*v*v*v*w*w*w*w;
				*(KURTuv+k) += u*u*u*u*v*v*v*v;
	
				*(dd+k) += dfl*dfl;
								
				if (j==0)
				{
					jp=yd-1;
					jpp=yd-2;
					jn=1;
					jnn=2;
				}
				else if (j==1)
				{
					jp=0;
					jpp=yd-1;
					jn=2;
					jnn=3;
				}
				else if (j==(yd-2))
				{
					jp=j-1;
					jpp=j-2;
					jn=j+1;
					jnn=0;
				}
				else if (j==(yd-1))
				{
					jp=j-1;
					jpp=j-2;
					jn=0;
					jnn=1;
				}
				else
				{
					jp=j-1;
					jpp=j-2;
					jn=j+1;
					jnn=j+2;
				}
				ip=i-1;
				ipp=i-2;
				in=i+1;
				inn=i+2;
				
				dudx=(2./3.)*(Fb(V.vel,in,j,k,0,xd,yd,zd,3)-Fb(V.vel,ip,j,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,0,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,0,xd,yd,zd,3));	
				dvdx=(2./3.)*(Fb(V.vel,in,j,k,1,xd,yd,zd,3)-Fb(V.vel,ip,j,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,1,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,1,xd,yd,zd,3));	
				dwdx=(2./3.)*(Fb(V.vel,in,j,k,2,xd,yd,zd,3)-Fb(V.vel,ip,j,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,inn,j,k,2,xd,yd,zd,3)-Fb(V.vel,ipp,j,k,2,xd,yd,zd,3));
				
				dudy=(2./3.)*(Fb(V.vel,i,jn,k,0,xd,yd,zd,3)-Fb(V.vel,i,jp,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,0,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,0,xd,yd,zd,3));	
				dvdy=(2./3.)*(Fb(V.vel,i,jn,k,1,xd,yd,zd,3)-Fb(V.vel,i,jp,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,1,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,1,xd,yd,zd,3));	
				dwdy=(2./3.)*(Fb(V.vel,i,jn,k,2,xd,yd,zd,3)-Fb(V.vel,i,jp,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,jnn,k,2,xd,yd,zd,3)-Fb(V.vel,i,jpp,k,2,xd,yd,zd,3));
				
				if (k==1)	
				{
					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),2,xd,yd,zd,3)-Fb(V.vel,i,j,k,2,xd,yd,zd,3));
					dwdz=-(7./3.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)+3.*tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),1,xd,yd,zd,3)-Fb(V.vel,i,j,k,1,xd,yd,zd,3));
					dvdz=-(7./3.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)+3.*tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+4),0,xd,yd,zd,3)-Fb(V.vel,i,j,k,0,xd,yd,zd,3));
					dudz=-(7./3.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)+6.*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-3.*Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-(2./3.)*Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)+3.*tmp;
				}
				else if (k==(zd-2))
				{
					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),2,xd,yd,zd,3));
					dwdz=(7./3.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3)+3.*tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),1,xd,yd,zd,3));
					dvdz=(7./3.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3)+3.*tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,k,0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-4),0,xd,yd,zd,3));
					dudz=(7./3.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)-6.*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)+3.*Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3)+(2./3.)*Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3)+3.*tmp;
				}
				else if (k==2)
				{
					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-Fb(V.vel,i,j,k,2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3));
					dwdz=-(1./6.)*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-Fb(V.vel,i,j,k,1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3));
					dvdz=-(1./6.)*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-Fb(V.vel,i,j,k,0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+3),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3));
					dudz=-(1./6.)*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)+(1./6.)*Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-tmp;
				}
				else if (k==(zd-3))
				{
					tmp=(2./3.)*(Fb(V.vel,i,j,k,2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),2,xd,yd,zd,3));
					dwdz=(1./6.)*Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,2,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3)-tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,k,1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),1,xd,yd,zd,3));
					dvdz=(1./6.)*Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,1,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3)-tmp;
					tmp=(2./3.)*(Fb(V.vel,i,j,k,0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-3),0,xd,yd,zd,3));
					dudz=(1./6.)*Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)+(3./2.)*Fb(V.vel,i,j,k,0,xd,yd,zd,3)-(3./2.)*Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3)-(1./6.)*Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3)-tmp;
				}
				else
				{
					dwdz=(2./3.)*(Fb(V.vel,i,j,(k+1),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),2,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),2,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),2,xd,yd,zd,3));
					dvdz=(2./3.)*(Fb(V.vel,i,j,(k+1),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),1,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),1,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),1,xd,yd,zd,3));
					dudz=(2./3.)*(Fb(V.vel,i,j,(k+1),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-1),0,xd,yd,zd,3))-(1./12.)*(Fb(V.vel,i,j,(k+2),0,xd,yd,zd,3)-Fb(V.vel,i,j,(k-2),0,xd,yd,zd,3));
				}
				
				KE2d[j][k] += 0.5*(u*u+v*v+w*w);
				uu2d[j][k] += u*u;
				uv2d[j][k] += u*v;
				uw2d[j][k] += u*w;
				vv2d[j][k] += v*v;
				vw2d[j][k] += v*w;
				ww2d[j][k] += w*w;
				
				usuv2d[j][k] += u*0.5*(dudy+dvdx);
				usuw2d[j][k] += u*0.5*(dudz+dwdx);
				usvv2d[j][k] += u*0.5*(dvdy+dvdy);
				usww2d[j][k] += u*0.5*(dwdz+dwdz);
				usvw2d[j][k] += u*0.5*(dvdz+dwdy);
				
				vsvv2d[j][k] += v*0.5*(dvdy+dvdy);
				vsuv2d[j][k] += v*0.5*(dudy+dvdx);
				vsuw2d[j][k] += v*0.5*(dudz+dwdx);
				vsvw2d[j][k] += v*0.5*(dvdz+dwdy);
				vsww2d[j][k] += v*0.5*(dwdz+dwdz);
				
				wsvv2d[j][k] += w*0.5*(dvdy+dvdy);
				wsww2d[j][k] += w*0.5*(dwdz+dwdz);
				wsuv2d[j][k] += w*0.5*(dudy+dvdx);
				wsuw2d[j][k] += w*0.5*(dudz+dwdx);
				wswv2d[j][k] += w*0.5*(dwdy+dvdz);
				
				vke2d[j][k] += v*0.5*(u*u+v*v+w*w);
				wke2d[j][k] += w*0.5*(u*u+v*v+w*w);
	
				pdudx2d[j][k] += dfl*dudx;
				pdudy2d[j][k] += dfl*dudy;
				pdudz2d[j][k] += dfl*dudz;
				
				pdvdx2d[j][k] += dfl*dvdx;
				pdvdy2d[j][k] += dfl*dvdy;
				pdvdz2d[j][k] += dfl*dvdz;
				
				pdwdz2d[j][k] += dfl*dwdz;
				pdwdy2d[j][k] += dfl*dwdy;
				pdwdx2d[j][k] += dfl*dwdx;

				uuv2d[j][k] += u*u*v;
				uuw2d[j][k] += u*u*w;
				vvv2d[j][k] += v*v*v;
				vvw2d[j][k] += v*v*w;
				wwv2d[j][k] += w*w*v;
				www2d[j][k] += w*w*w;
				uvv2d[j][k] += u*v*v;
				uvw2d[j][k] += u*v*w;
				uwv2d[j][k] += u*v*w;
				uww2d[j][k] += u*w*w;
				vwv2d[j][k] += v*w*v;
				vww2d[j][k] += v*w*w;
				
				pu2d[j][k] += dfl*u;
				pv2d[j][k] += dfl*v;
				pw2d[j][k] += dfl*w;
				
				Epuu2d[j][k] -= 2.*nu*(dudx*(dudx+dudx)+dudy*(dudy+dvdx)+dudz*(dudz+dwdx));
				Epvv2d[j][k] -= 2.*nu*(dvdx*(dvdx+dudy)+dvdy*(dvdy+dvdy)+dvdz*(dvdz+dwdy));
				Epww2d[j][k] -= 2.*nu*(dwdx*(dwdx+dudz)+dwdy*(dwdy+dvdz)+dwdz*(dwdz+dwdz));
				Epuw2d[j][k] -= nu*(dudx*(dwdx+dudz)+dudy*(dwdy+dvdz)+dudz*(dwdz+dwdz)+dwdx*(dudx+dudx)+dwdy*(dudy+dvdx)+dwdz*(dudz+dwdx));
				Epuv2d[j][k] -= nu*(dudx*(dudy+dvdx)+dudy*(dvdy+dvdy)+dudz*(dvdz+dwdy)+dvdx*(dudx+dudx)+dvdy*(dudy+dvdx)+dvdz*(dudz+dwdx));
				Epvw2d[j][k] -= nu*(dvdx*(dwdx+dudz)+dvdy*(dvdz+dwdy)+dvdz*(dwdz+dwdz)+dwdx*(dvdx+dudy)+dwdy*(dvdy+dvdy)+dwdz*(dvdz+dwdy));
				Epii2d[j][k] -= nu*((dudx*(dudx+dudx)+dudy*(dudy+dvdx)+dudz*(dudz+dwdx))+(dvdx*(dvdx+dudy)+dvdy*(dvdy+dvdy)+dvdz*(dvdz+dwdy))+(dwdx*(dwdx+dudz)+dwdy*(dwdy+dvdz)+dwdz*(dwdz+dwdz)));
				
				ox=fac*(dwdy-dvdz);
				oy=fac*(dudz-dwdx);
				oz=fac*(dvdx-dudy);

				*(wx+k) += ox*ox;
             	*(wy+k) += oy*oy;
           		*(wz+k) += oz*oz;
				
				omegax2d[j][k] += ox*ox;
				omegay2d[j][k] += oy*oy;
				omegaz2d[j][k] += oz*oz;

				dd2d[j][k] += dfl*dfl;
				
				if (u>0. && w>0.)
				{
				  quadrant[0][k] += u*w;
				  nquadrant[0][k] ++;
				  
				  quadrant2d[0][j][k] += u*w;
				}
				else if (u<0. && w>0.)
				{
				  quadrant[1][k] += u*w;
				  nquadrant[1][k] ++;
				  
				  quadrant2d[1][j][k] += u*w;
				}
				else if (u<0. && w<0.)
				{
				  quadrant[2][k] += u*w;
				  nquadrant[2][k] ++;
				  
				  quadrant2d[2][j][k] += u*w;
				}
				else if (u>0. && w<0.)
				{
				  quadrant[3][k] += u*w;
				  nquadrant[3][k] ++;
				  
				  quadrant2d[3][j][k] += u*w;
				}
			}
		}
	}
	
	for (j=0;j<exty;j++)
	{
		if (j==0)
		{
			jp=yd-1;
			jpp=yd-2;
			jn=1;
			jnn=2;
		}
		else if (j==1)
		{
			jp=0;
			jpp=yd-1;
			jn=2;
			jnn=3;
		}
		else if (j==(yd-2))
		{
			jp=j-1;
			jpp=j-2;
			jn=j+1;
			jnn=0;
		}
		else if (j==(yd-1))
		{
			jp=j-1;
			jpp=j-2;
			jn=0;
			jnn=1;
		}
		else
		{
			jp=j-1;
			jpp=j-2;
			jn=j+1;
			jnn=j+2;
		}
		for (k=1;k<(zd-1);k++)
		{
			
			dUdy2d[j][k] = (2./3.)*(uAVG[jn][k]-uAVG[jp][k])-(1./12.)*(uAVG[jnn][k]-uAVG[jpp][k]);
			dVdy2d[j][k] = (2./3.)*(vAVG[jn][k]-vAVG[jp][k])-(1./12.)*(vAVG[jnn][k]-vAVG[jpp][k]);
			dWdy2d[j][k] = (2./3.)*(wAVG[jn][k]-wAVG[jp][k])-(1./12.)*(wAVG[jnn][k]-wAVG[jpp][k]);
			
			if (k==1)	
			{
				tmp = (2./3.)*(uAVG[j][k+3]-uAVG[j][k+1])-(1./12.)*(uAVG[j][k+4]-uAVG[j][k]);
				dUdz2d[j][k] = -(7./3.)*uAVG[j][k]+6.*uAVG[j][k+1]-3.*uAVG[j][k+2]-(2./3.)*uAVG[j][k+3]+3.*tmp;
				tmp = (2./3.)*(vAVG[j][k+3]-vAVG[j][k+1])-(1./12.)*(vAVG[j][k+4]-vAVG[j][k]);
				dVdz2d[j][k] = -(7./3.)*vAVG[j][k]+6.*vAVG[j][k+1]-3.*vAVG[j][k+2]-(2./3.)*vAVG[j][k+3]+3.*tmp;
				tmp = (2./3.)*(wAVG[j][k+3]-wAVG[j][k+1])-(1./12.)*(wAVG[j][k+4]-wAVG[j][k]);
				dWdz2d[j][k] = -(7./3.)*wAVG[j][k]+6.*wAVG[j][k+1]-3.*wAVG[j][k+2]-(2./3.)*wAVG[j][k+3]+3.*tmp;
			}
			else if (k==(zd-2))
			{
				tmp = (2./3.)*(uAVG[j][k-1]-uAVG[j][k-3])-(1./12.)*(uAVG[j][k]-uAVG[j][k-4]);
				dUdz2d[j][k] = (7./3.)*uAVG[j][k]-6.*uAVG[j][k-1]+3.*uAVG[j][k-2]+(2./3.)*uAVG[j][k-3]+3.*tmp;
				tmp = (2./3.)*(vAVG[j][k-1]-vAVG[j][k-3])-(1./12.)*(vAVG[j][k]-vAVG[j][k-4]);
				dVdz2d[j][k] = (7./3.)*vAVG[j][k]-6.*vAVG[j][k-1]+3.*vAVG[j][k-2]+(2./3.)*vAVG[j][k-3]+3.*tmp;
				tmp = (2./3.)*(wAVG[j][k-1]-wAVG[j][k-3])-(1./12.)*(wAVG[j][k]-wAVG[j][k-4]);
				dWdz2d[j][k] = (7./3.)*wAVG[j][k]-6.*wAVG[j][k-1]+3.*wAVG[j][k-2]+(2./3.)*wAVG[j][k-3]+3.*tmp;
			}
			else if (k==2)
			{
				tmp = (2./3.)*(uAVG[j][k+2]-uAVG[j][k])-(1./12.)*(uAVG[j][k+3]-uAVG[j][k-1]);
				dUdz2d[j][k] = -(1./6.)*uAVG[j][k-1]-(3./2.)*uAVG[j][k]+(3./2.)*uAVG[j][k+1]+(1./6.)*uAVG[j][k+2]-tmp;
				tmp = (2./3.)*(vAVG[j][k+2]-vAVG[j][k])-(1./12.)*(vAVG[j][k+3]-vAVG[j][k-1]);
				dVdz2d[j][k] = -(1./6.)*vAVG[j][k-1]-(3./2.)*vAVG[j][k]+(3./2.)*vAVG[j][k+1]+(1./6.)*vAVG[j][k+2]-tmp;
				tmp = (2./3.)*(wAVG[j][k+2]-wAVG[j][k])-(1./12.)*(wAVG[j][k+3]-wAVG[j][k-1]);
				dWdz2d[j][k] = -(1./6.)*wAVG[j][k-1]-(3./2.)*wAVG[j][k]+(3./2.)*wAVG[j][k+1]+(1./6.)*wAVG[j][k+2]-tmp;
			}
			else if (k==(zd-3))
			{
				tmp = (2./3.)*(uAVG[j][k]-uAVG[j][k-2])-(1./12.)*(uAVG[j][k+1]-uAVG[j][k-3]);
				dUdz2d[j][k] = (1./6.)*uAVG[j][k+1]+(3./2.)*uAVG[j][k]-(3./2.)*uAVG[j][k-1]-(1./6.)*uAVG[j][k-2]-tmp;
				tmp = (2./3.)*(vAVG[j][k]-vAVG[j][k-2])-(1./12.)*(vAVG[j][k+1]-vAVG[j][k-3]);
				dVdz2d[j][k] = (1./6.)*vAVG[j][k+1]+(3./2.)*vAVG[j][k]-(3./2.)*vAVG[j][k-1]-(1./6.)*vAVG[j][k-2]-tmp;
				tmp = (2./3.)*(wAVG[j][k]-wAVG[j][k-2])-(1./12.)*(wAVG[j][k+1]-wAVG[j][k-3]);
				dWdz2d[j][k] = (1./6.)*wAVG[j][k+1]+(3./2.)*wAVG[j][k]-(3./2.)*wAVG[j][k-1]-(1./6.)*wAVG[j][k-2]-tmp;
			}
			else
			{
				dUdz2d[j][k] = (2./3.)*(uAVG[j][k+1]-uAVG[j][k-1])-(1./12.)*(uAVG[j][k+2]-uAVG[j][k-2]);
				dVdz2d[j][k] = (2./3.)*(vAVG[j][k+1]-vAVG[j][k-1])-(1./12.)*(vAVG[j][k+2]-vAVG[j][k-2]);
				dWdz2d[j][k] = (2./3.)*(wAVG[j][k+1]-wAVG[j][k-1])-(1./12.)*(wAVG[j][k+2]-wAVG[j][k-2]);
			}
		}
	}
	
	for (j=0;j<exty;j++)
	{
		for (k=0;k<zd;k++)
		{
			OmegaX[j][k] = fac*(dWdy2d[j][k]-dVdz2d[j][k]);
		}
	}
	
	for (j=0;j<exty;j++)
	{
		for (k=1;k<(zd-1);k++)
		{
			KE2d[j][k] /= ((double)(pd.numproc*extx));
			uu2d[j][k] /= ((double)(pd.numproc*extx));
			uv2d[j][k] /= ((double)(pd.numproc*extx));
			uw2d[j][k] /= ((double)(pd.numproc*extx));
			vv2d[j][k] /= ((double)(pd.numproc*extx));
			vw2d[j][k] /= ((double)(pd.numproc*extx));
			ww2d[j][k] /= ((double)(pd.numproc*extx));
			
			usuv2d[j][k] /= ((double)(pd.numproc*extx));
			usuw2d[j][k] /= ((double)(pd.numproc*extx));
			vsvv2d[j][k] /= ((double)(pd.numproc*extx));
			wswv2d[j][k] /= ((double)(pd.numproc*extx));
			usvw2d[j][k] /= ((double)(pd.numproc*extx));
			vsvw2d[j][k] /= ((double)(pd.numproc*extx));
			vsww2d[j][k] /= ((double)(pd.numproc*extx));
			wsww2d[j][k] /= ((double)(pd.numproc*extx));
			wsuv2d[j][k] /= ((double)(pd.numproc*extx));
			usww2d[j][k] /= ((double)(pd.numproc*extx));
			vsuv2d[j][k] /= ((double)(pd.numproc*extx));
			wsvv2d[j][k] /= ((double)(pd.numproc*extx));
			wsuw2d[j][k] /= ((double)(pd.numproc*extx));
			usvv2d[j][k] /= ((double)(pd.numproc*extx));
			vsuw2d[j][k] /= ((double)(pd.numproc*extx));
			
			pdudx2d[j][k] /= ((double)(pd.numproc*extx));
			pdvdy2d[j][k] /= ((double)(pd.numproc*extx));
			pdwdz2d[j][k] /= ((double)(pd.numproc*extx));
			pdwdy2d[j][k] /= ((double)(pd.numproc*extx));
			pdvdx2d[j][k] /= ((double)(pd.numproc*extx));
			pdudy2d[j][k] /= ((double)(pd.numproc*extx));
			pdudz2d[j][k] /= ((double)(pd.numproc*extx));
			pdwdx2d[j][k] /= ((double)(pd.numproc*extx));
			pdvdz2d[j][k] /= ((double)(pd.numproc*extx));
				
			vke2d[j][k] /= ((double)(pd.numproc*extx));
			wke2d[j][k] /= ((double)(pd.numproc*extx));
			
			uuv2d[j][k] /= ((double)(pd.numproc*extx));
			uuw2d[j][k] /= ((double)(pd.numproc*extx));
			vvv2d[j][k] /= ((double)(pd.numproc*extx));
			vvw2d[j][k] /= ((double)(pd.numproc*extx));
			wwv2d[j][k] /= ((double)(pd.numproc*extx));
			www2d[j][k] /= ((double)(pd.numproc*extx));
			uvv2d[j][k] /= ((double)(pd.numproc*extx));
			uvw2d[j][k] /= ((double)(pd.numproc*extx));
			uwv2d[j][k] /= ((double)(pd.numproc*extx));
			uww2d[j][k] /= ((double)(pd.numproc*extx));
			vwv2d[j][k] /= ((double)(pd.numproc*extx));
			vww2d[j][k] /= ((double)(pd.numproc*extx));
				
			Epii2d[j][k] /= ((double)(pd.numproc*extx));
			Epuu2d[j][k] /= ((double)(pd.numproc*extx));
			Epvv2d[j][k] /= ((double)(pd.numproc*extx));
			Epww2d[j][k] /= ((double)(pd.numproc*extx));
			Epuw2d[j][k] /= ((double)(pd.numproc*extx));
			Epuv2d[j][k] /= ((double)(pd.numproc*extx));
			Epvw2d[j][k] /= ((double)(pd.numproc*extx));
			
			pu2d[j][k] /= ((double)(pd.numproc*extx));
			pv2d[j][k] /= ((double)(pd.numproc*extx));
			pw2d[j][k] /= ((double)(pd.numproc*extx));
			
			omegax2d[j][k] /= ((double)(pd.numproc*extx));
			omegay2d[j][k] /= ((double)(pd.numproc*extx));
			omegaz2d[j][k] /= ((double)(pd.numproc*extx));

			dd2d[j][k] /= ((double)(pd.numproc*extx));
			
			quadrant2d[0][j][k] /= ((double)(pd.numproc*extx));
			quadrant2d[1][j][k] /= ((double)(pd.numproc*extx));
			quadrant2d[2][j][k] /= ((double)(pd.numproc*extx));
			quadrant2d[3][j][k] /= ((double)(pd.numproc*extx));
		}
	}
	
	for (k=0;k<zd;k++)
	{
		*(uu+k) /= (double)(num);
		*(uv+k) /= (double)(num);
		*(uw+k) /= (double)(num);
		*(vv+k) /= (double)(num);
		*(vw+k) /= (double)(num);
		*(ww+k) /= (double)(num);
		*(dd+k) /= (double)(num);
		
		*(SKEWu+k) /= (double)(num);
		*(SKEWv+k) /= (double)(num);
		*(SKEWw+k) /= (double)(num);
		*(SKEWp+k) /= (double)(num);
		*(SKEWuw+k) /= (double)(num);
		*(SKEWuv+k) /= (double)(num);
		*(SKEWvw+k) /= (double)(num);
		
		*(KURTu+k) /= (double)(num);
		*(KURTv+k) /= (double)(num);
		*(KURTw+k) /= (double)(num);
		*(KURTp+k) /= (double)(num);
		*(KURTuw+k) /= (double)(num);
		*(KURTuv+k) /= (double)(num);
		*(KURTvw+k) /= (double)(num);

		*(wx+k) /= (double)(num);
     		*(wy+k) /= (double)(num);
     		*(wz+k) /= (double)(num);
		
		quadrant[0][k] /= (double)(num);
		quadrant[1][k] /= (double)(num);
		quadrant[2][k] /= (double)(num);
		quadrant[3][k] /= (double)(num);
	}
	
	MPI_Reduce(uu,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uu,DUMMY,zd*sizeof(DP));

	MPI_Reduce(uv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uv,DUMMY,zd*sizeof(DP));

	MPI_Reduce(uw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uw,DUMMY,zd*sizeof(DP));

	MPI_Reduce(vv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vv,DUMMY,zd*sizeof(DP));

	MPI_Reduce(vw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vw,DUMMY,zd*sizeof(DP));

	MPI_Reduce(ww,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(ww,DUMMY,zd*sizeof(DP));

	MPI_Reduce(dd,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(dd,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWu,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWu,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWv,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWp,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWp,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWuw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWuw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWuv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWuv,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(SKEWvw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(SKEWvw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTu,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTu,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTv,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTp,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTp,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTuw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTuw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTvw,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTvw,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(KURTuv,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KURTuv,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(wx,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wx,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(wy,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wy,DUMMY,zd*sizeof(DP));
	
	MPI_Reduce(wz,DUMMY,zd,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMY,zd,MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wz,DUMMY,zd*sizeof(DP));
	
	MPI_Allreduce(&quadrant[0][0],DUMMYq,4*zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(&quadrant[0][0],DUMMYq,4*zd*sizeof(DP));
	
	MPI_Allreduce(&nquadrant[0][0],iDUMMYq,4*zd,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
	memcpy(&nquadrant[0][0],iDUMMYq,4*zd*sizeof(int));
	
/*****************************************************************************************
 * 2D part                                                                               *
 *****************************************************************************************/
 	
	MPI_Reduce(&KE2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(KE2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&uu2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uu2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&uv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&uw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&vv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&vw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&ww2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(ww2d,DUMMYYZ,yd*zd*sizeof(DP));
	
//************************************************************************************************************
	
	MPI_Reduce(&usuv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(usuv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&usuw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(usuw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&usvv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(usvv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&usww2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(usww2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&usvw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(usvw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	
	MPI_Reduce(&vsvv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vsvv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&vsuv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vsuv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&vsuw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vsuw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&vsvw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vsvw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&vsww2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vsww2d,DUMMYYZ,yd*zd*sizeof(DP));

	
	MPI_Reduce(&wsvv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wsvv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&wsww2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wsww2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&wsuv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wsuv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&wsuw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wsuw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&wswv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wswv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
//*************************************************************************************************************
	
	MPI_Reduce(&vke2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vke2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&wke2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wke2d,DUMMYYZ,yd*zd*sizeof(DP));
	
//*************************************************************************************************************
	
	MPI_Reduce(&Epii2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epii2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&Epuu2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epuu2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&Epvv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epvv2d,DUMMYYZ,yd*zd*sizeof(DP));
		
	MPI_Reduce(&Epww2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epww2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&Epuv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epuv2d,DUMMYYZ,yd*zd*sizeof(DP));	
		
	MPI_Reduce(&Epuw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epuw2d,DUMMYYZ,yd*zd*sizeof(DP));

	MPI_Reduce(&Epvw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(Epvw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
//**************************************************************************************************************
	
	MPI_Reduce(&pdudx2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pdudx2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&pdudy2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pdudy2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&pdudz2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pdudz2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&pdvdx2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pdvdx2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&pdvdy2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pdvdy2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&pdvdz2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pdvdz2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&pdwdx2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pdwdx2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&pdwdy2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pdwdy2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&pdwdz2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pdwdz2d,DUMMYYZ,yd*zd*sizeof(DP));
	
//***************************************************************************************************************
	
	MPI_Reduce(&uuv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uuv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&uuw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uuw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&vvv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vvv2d,DUMMYYZ,yd*zd*sizeof(DP));

	MPI_Reduce(&vvw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vvw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&wwv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(wwv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&www2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(www2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&uvv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uvv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&uvw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uvw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&uwv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uwv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&uww2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(uww2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&vwv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vwv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&vww2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(vww2d,DUMMYYZ,yd*zd*sizeof(DP));
	
//***************************************************************************************************************
	
	MPI_Reduce(&pu2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pu2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&pv2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pv2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&pw2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(pw2d,DUMMYYZ,yd*zd*sizeof(DP));
	
//****************************************************************************************************************

	MPI_Reduce(&omegax2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(omegax2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&omegay2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(omegay2d,DUMMYYZ,yd*zd*sizeof(DP));
	
	MPI_Reduce(&omegaz2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
	memcpy(omegaz2d,DUMMYYZ,yd*zd*sizeof(DP));
	
//****************************************************************************************************************

	MPI_Reduce(&dd2d[0][0],DUMMYYZ,(yd*zd),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Bcast(DUMMYYZ,(yd*zd),MPI_DOUBLE,0,MPI_COMM_WORLD);
    memcpy(dd2d,DUMMYYZ,yd*zd*sizeof(DP));	
	
	MPI_Allreduce(&quadrant2d[0][0][0],DUMMYq2d,4*yd*zd,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	memcpy(&quadrant2d[0][0][0],DUMMYq2d,4*yd*zd*sizeof(DP));	

//****************************************************************************************************************

	for (j=0;j<exty;j++)
	{
		if (j==0)
		{
			jp=yd-1;
			jpp=yd-2;
			jn=1;
			jnn=2;
		}
		else if (j==1)
		{
			jp=0;
			jpp=yd-1;
			jn=2;
			jnn=3;
		}
		else if (j==(yd-2))
		{
			jp=j-1;
			jpp=j-2;
			jn=j+1;
			jnn=0;
		}
		else if (j==(yd-1))
		{
			jp=j-1;
			jpp=j-2;
			jn=0;
			jnn=1;
		}
		else
		{
			jp=j-1;
			jpp=j-2;
			jn=j+1;
			jnn=j+2;
		}
		
		for (k=1;k<(zd-1);k++)
		{
			Vduudy2d[j][k] = vAVG[j][k]*((2./3.)*(uu2d[jn][k]-uu2d[jp][k])-(1./12.)*(uu2d[jnn][k]-uu2d[jpp][k]));
			Vdvvdy2d[j][k] = vAVG[j][k]*((2./3.)*(vv2d[jn][k]-vv2d[jp][k])-(1./12.)*(vv2d[jnn][k]-vv2d[jpp][k]));
			Vdwwdy2d[j][k] = vAVG[j][k]*((2./3.)*(ww2d[jn][k]-ww2d[jp][k])-(1./12.)*(ww2d[jnn][k]-ww2d[jpp][k]));
			Vduvdy2d[j][k] = vAVG[j][k]*((2./3.)*(uv2d[jn][k]-uv2d[jp][k])-(1./12.)*(uv2d[jnn][k]-uv2d[jpp][k]));
			Vduwdy2d[j][k] = vAVG[j][k]*((2./3.)*(uw2d[jn][k]-uw2d[jp][k])-(1./12.)*(uw2d[jnn][k]-uw2d[jpp][k]));
			Vdvwdy2d[j][k] = vAVG[j][k]*((2./3.)*(vw2d[jn][k]-vw2d[jp][k])-(1./12.)*(vw2d[jnn][k]-vw2d[jpp][k]));

			VdKdy2d[j][k] = vAVG[j][k]*((2./3.)*(KE2d[jn][k]-KE2d[jp][k])-(1./12.)*(KE2d[jnn][k]-KE2d[jpp][k]));
			
			Pii2d[j][k] = -(uv2d[j][k]*dUdy2d[j][k]+uw2d[j][k]*dUdz2d[j][k]+vv2d[j][k]*dVdy2d[j][k]+vw2d[j][k]*dVdz2d[j][k]+vw2d[j][k]*dWdy2d[j][k]+ww2d[j][k]*dWdz2d[j][k]);
			P1[j][k] = - uv2d[j][k]*dUdy2d[j][k];
			P2[j][k] = - uw2d[j][k]*dUdz2d[j][k];
			P3[j][k] = - vv2d[j][k]*dVdy2d[j][k];
			P4[j][k] = - vw2d[j][k]*dVdz2d[j][k];
			P5[j][k] = - vw2d[j][k]*dWdy2d[j][k];
			P6[j][k] = - ww2d[j][k]*dWdz2d[j][k];
			P7[j][k] = - uv2d[j][k]*dWdy2d[j][k];
			P8[j][k] = - uw2d[j][k]*dWdz2d[j][k];
			P9[j][k] = - vw2d[j][k]*dUdy2d[j][k];
			P10[j][k] = - ww2d[j][k]*dUdz2d[j][k];
			P11[j][k] = - uv2d[j][k]*dVdy2d[j][k];
			P12[j][k] = - uw2d[j][k]*dVdz2d[j][k];
			P13[j][k] = - vv2d[j][k]*dUdy2d[j][k];
			P14[j][k] = - vw2d[j][k]*dUdz2d[j][k];
			P15[j][k] = - vv2d[j][k]*dWdy2d[j][k];
			P16[j][k] = - vw2d[j][k]*dWdz2d[j][k];
			P17[j][k] = - vw2d[j][k]*dVdy2d[j][k];
			P18[j][k] = - ww2d[j][k]*dVdz2d[j][k];
			Puu2d[j][k] = -2.*(uv2d[j][k]*dUdy2d[j][k]+uw2d[j][k]*dUdz2d[j][k]);
			Pvv2d[j][k] = -2.*(vv2d[j][k]*dVdy2d[j][k]+vw2d[j][k]*dVdz2d[j][k]);
			Pww2d[j][k] = -2.*(vw2d[j][k]*dWdy2d[j][k]+ww2d[j][k]*dWdz2d[j][k]);
			Puw2d[j][k] = -(uv2d[j][k]*dWdy2d[j][k]+uw2d[j][k]*dWdz2d[j][k]+vw2d[j][k]*dUdy2d[j][k]+ww2d[j][k]*dUdz2d[j][k]);
			Puv2d[j][k] = -(uv2d[j][k]*dVdy2d[j][k]+uw2d[j][k]*dVdz2d[j][k]+vv2d[j][k]*dUdy2d[j][k]+vw2d[j][k]*dUdz2d[j][k]);
			Pvw2d[j][k] = -(vv2d[j][k]*dWdy2d[j][k]+vw2d[j][k]*dWdz2d[j][k]+vw2d[j][k]*dVdy2d[j][k]+ww2d[j][k]*dVdz2d[j][k]);
			
			Paiuu2d[j][k] = 2.*pdudx2d[j][k];
			Paivv2d[j][k] = 2.*pdvdy2d[j][k];
			Paiww2d[j][k] = 2.*pdwdz2d[j][k];
			Paiuv2d[j][k] = pdudy2d[j][k]+pdvdx2d[j][k];
			Paiuw2d[j][k] = pdudz2d[j][k]+pdwdx2d[j][k];
			Paivw2d[j][k] = pdwdy2d[j][k]+pdvdz2d[j][k];
			
			Tiipai2d[j][k] = -0.5*2.*((2./3.)*(pv2d[jn][k]-pv2d[jp][k])-(1./12.)*(pv2d[jnn][k]-pv2d[jpp][k]));
			Tuupai2d[j][k] = 0.;
			Tvvpai2d[j][k] = -2.*((2./3.)*(pv2d[jn][k]-pv2d[jp][k])-(1./12.)*(pv2d[jnn][k]-pv2d[jpp][k]));
			Twwpai2d[j][k] = 0.;
			Tuvpai2d[j][k] = -((2./3.)*(pu2d[jn][k]-pu2d[jp][k])-(1./12.)*(pu2d[jnn][k]-pu2d[jpp][k]));
			Tuwpai2d[j][k] = 0.;
			Tvwpai2d[j][k] = -((2./3.)*(pw2d[jn][k]-pw2d[jp][k])-(1./12.)*(pw2d[jnn][k]-pw2d[jpp][k]));
			
			Tiis2d[j][k] = 2.*nu*((2./3.)*(usuv2d[jn][k]-usuv2d[jp][k]+vsvv2d[jn][k]-vsvv2d[jp][k]+wswv2d[jn][k]-wswv2d[jp][k])-(1./12.)*(usuv2d[jnn][k]-usuv2d[jpp][k]+vsvv2d[jnn][k]-vsvv2d[jpp][k]+wswv2d[jnn][k]-wswv2d[jpp][k]));
			Tuus2d[j][k] = 4.*nu*((2./3.)*(usuv2d[jn][k]-usuv2d[jp][k])-(1./12.)*(usuv2d[jnn][k]-usuv2d[jpp][k]));
			Tvvs2d[j][k] = 4.*nu*((2./3.)*(vsvv2d[jn][k]-vsvv2d[jp][k])-(1./12.)*(vsvv2d[jnn][k]-vsvv2d[jpp][k]));
			Twws2d[j][k] = 4.*nu*((2./3.)*(wswv2d[jn][k]-wswv2d[jp][k])-(1./12.)*(wswv2d[jnn][k]-wswv2d[jpp][k]));
			Tuvs2d[j][k] = 2.*nu*((2./3.)*(vsuv2d[jn][k]-vsuv2d[jp][k]+usvv2d[jn][k]-usvv2d[jp][k])-(1./12.)*(vsuv2d[jnn][k]-vsuv2d[jpp][k]+usvv2d[jnn][k]-usvv2d[jpp][k]));
			Tuws2d[j][k] = 2.*nu*((2./3.)*(wsuv2d[jn][k]-wsuv2d[jp][k]+usvw2d[jn][k]-usvw2d[jp][k])-(1./12.)*(wsuv2d[jnn][k]-wsuv2d[jpp][k]+usvw2d[jnn][k]-usvw2d[jpp][k]));
			Tvws2d[j][k] = 2.*nu*((2./3.)*(wsvv2d[jn][k]-wsvv2d[jp][k]+vsvw2d[jn][k]-vsvw2d[jp][k])-(1./12.)*(wsvv2d[jnn][k]-wsvv2d[jpp][k]+vsvw2d[jnn][k]-vsvw2d[jpp][k]));
			
			Tiir2d[j][k] = -((2./3.)*(vke2d[jn][k]-vke2d[jp][k])-(1./12.)*(vke2d[jnn][k]-vke2d[jpp][k]));
			Tuur2d[j][k] = -((2./3.)*(uuv2d[jn][k]-uuv2d[jp][k])-(1./12.)*(uuv2d[jnn][k]-uuv2d[jpp][k]));
			Tvvr2d[j][k] = -((2./3.)*(vvv2d[jn][k]-vvv2d[jp][k])-(1./12.)*(vvv2d[jnn][k]-vvv2d[jpp][k]));
			Twwr2d[j][k] = -((2./3.)*(wwv2d[jn][k]-wwv2d[jp][k])-(1./12.)*(wwv2d[jnn][k]-wwv2d[jpp][k]));
			Tuvr2d[j][k] = -((2./3.)*(uvv2d[jn][k]-uvv2d[jp][k])-(1./12.)*(uvv2d[jnn][k]-uvv2d[jpp][k]));
			Tuwr2d[j][k] = -((2./3.)*(uwv2d[jn][k]-uwv2d[jp][k])-(1./12.)*(uwv2d[jnn][k]-uwv2d[jpp][k]));
			Tvwr2d[j][k] = -((2./3.)*(vwv2d[jn][k]-vwv2d[jp][k])-(1./12.)*(vwv2d[jnn][k]-vwv2d[jpp][k]));
			
			if (k==1)	
			{
				tmp = (2./3.)*(uu2d[j][k+3]-uu2d[j][k+1])-(1./12.)*(uu2d[j][k+4]-uu2d[j][k]);
				Wduudz2d[j][k] = wAVG[j][k]*(-(7./3.)*uu2d[j][k]+6.*uu2d[j][k+1]-3.*uu2d[j][k+2]-(2./3.)*uu2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(vv2d[j][k+3]-vv2d[j][k+1])-(1./12.)*(vv2d[j][k+4]-vv2d[j][k]);
				Wdvvdz2d[j][k] = wAVG[j][k]*(-(7./3.)*vv2d[j][k]+6.*vv2d[j][k+1]-3.*vv2d[j][k+2]-(2./3.)*vv2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(ww2d[j][k+3]-ww2d[j][k+1])-(1./12.)*(ww2d[j][k+4]-ww2d[j][k]);
				Wdwwdz2d[j][k] = wAVG[j][k]*(-(7./3.)*ww2d[j][k]+6.*ww2d[j][k+1]-3.*ww2d[j][k+2]-(2./3.)*ww2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(uv2d[j][k+3]-uv2d[j][k+1])-(1./12.)*(uv2d[j][k+4]-uv2d[j][k]);
				Wduvdz2d[j][k] = wAVG[j][k]*(-(7./3.)*uv2d[j][k]+6.*uv2d[j][k+1]-3.*uv2d[j][k+2]-(2./3.)*uv2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(uw2d[j][k+3]-uw2d[j][k+1])-(1./12.)*(uw2d[j][k+4]-uw2d[j][k]);
				Wduwdz2d[j][k] = wAVG[j][k]*(-(7./3.)*uw2d[j][k]+6.*uw2d[j][k+1]-3.*uw2d[j][k+2]-(2./3.)*uw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(vw2d[j][k+3]-vw2d[j][k+1])-(1./12.)*(vw2d[j][k+4]-vw2d[j][k]);
				Wdvwdz2d[j][k] = wAVG[j][k]*(-(7./3.)*vw2d[j][k]+6.*vw2d[j][k+1]-3.*vw2d[j][k+2]-(2./3.)*vw2d[j][k+3]+3.*tmp);

				tmp = (2./3.)*(KE2d[j][k+3]-KE2d[j][k+1])-(1./12.)*(KE2d[j][k+4]-KE2d[j][k]);
				WdKdz2d[j][k] = wAVG[j][k]*(-(7./3.)*KE2d[j][k]+6.*KE2d[j][k+1]-3.*KE2d[j][k+2]-(2./3.)*KE2d[j][k+3]+3.*tmp);
				
				tmp = (2./3.)*(pw2d[j][k+3]-pw2d[j][k+1])-(1./12.)*(pw2d[j][k+4]-pw2d[j][k]);
				Tiipai2d[j][k] -= 0.5*2.*(-(7./3.)*pw2d[j][k]+6.*pw2d[j][k+1]-3.*pw2d[j][k+2]-(2./3.)*pw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(pw2d[j][k+3]-pw2d[j][k+1])-(1./12.)*(pw2d[j][k+4]-pw2d[j][k]);
				Twwpai2d[j][k] = -2.*(-(7./3.)*pw2d[j][k]+6.*pw2d[j][k+1]-3.*pw2d[j][k+2]-(2./3.)*pw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(pu2d[j][k+3]-pu2d[j][k+1])-(1./12.)*(pu2d[j][k+4]-pu2d[j][k]);
				Tuwpai2d[j][k] = -(-(7./3.)*pu2d[j][k]+6.*pu2d[j][k+1]-3.*pu2d[j][k+2]-(2./3.)*pu2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(pv2d[j][k+3]-pv2d[j][k+1])-(1./12.)*(pv2d[j][k+4]-pv2d[j][k]);
				Tvwpai2d[j][k] -= (-(7./3.)*pv2d[j][k]+6.*pv2d[j][k+1]-3.*pv2d[j][k+2]-(2./3.)*pv2d[j][k+3]+3.*tmp);
				
				tmp = (2./3.)*(usuw2d[j][k+3]-usuw2d[j][k+1])-(1./12.)*(usuw2d[j][k+4]-usuw2d[j][k]);
				Tiis2d[j][k] += 2.*nu*(-(7./3.)*usuw2d[j][k]+6.*usuw2d[j][k+1]-3.*usuw2d[j][k+2]-(2./3.)*usuw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(vsvw2d[j][k+3]-vsvw2d[j][k+1])-(1./12.)*(vsvw2d[j][k+4]-vsvw2d[j][k]);
				Tiis2d[j][k] += 2.*nu*(-(7./3.)*vsvw2d[j][k]+6.*vsvw2d[j][k+1]-3.*vsvw2d[j][k+2]-(2./3.)*vsvw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(wsww2d[j][k+3]-wsww2d[j][k+1])-(1./12.)*(wsww2d[j][k+4]-wsww2d[j][k]);
				Tiis2d[j][k] += 2.*nu*(-(7./3.)*wsww2d[j][k]+6.*wsww2d[j][k+1]-3.*wsww2d[j][k+2]-(2./3.)*wsww2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(usuw2d[j][k+3]-usuw2d[j][k+1])-(1./12.)*(usuw2d[j][k+4]-usuw2d[j][k]);
				Tuus2d[j][k] += 4.*nu*(-(7./3.)*usuw2d[j][k]+6.*usuw2d[j][k+1]-3.*usuw2d[j][k+2]-(2./3.)*usuw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(vsvw2d[j][k+3]-vsvw2d[j][k+1])-(1./12.)*(vsvw2d[j][k+4]-vsvw2d[j][k]);
				Tvvs2d[j][k] += 4.*nu*(-(7./3.)*vsvw2d[j][k]+6.*vsvw2d[j][k+1]-3.*vsvw2d[j][k+2]-(2./3.)*vsvw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(wsww2d[j][k+3]-wsww2d[j][k+1])-(1./12.)*(wsww2d[j][k+4]-wsww2d[j][k]);
				Twws2d[j][k] += 4.*nu*(-(7./3.)*wsww2d[j][k]+6.*wsww2d[j][k+1]-3.*wsww2d[j][k+2]-(2./3.)*wsww2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(vsuw2d[j][k+3]-vsuw2d[j][k+1])-(1./12.)*(vsuw2d[j][k+4]-vsuw2d[j][k]);
				Tuvs2d[j][k] += 2.*nu*(-(7./3.)*vsuw2d[j][k]+6.*vsuw2d[j][k+1]-3.*vsuw2d[j][k+2]-(2./3.)*vsuw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(usvw2d[j][k+3]-usvw2d[j][k+1])-(1./12.)*(usvw2d[j][k+4]-usvw2d[j][k]);
				Tuvs2d[j][k] += 2.*nu*(-(7./3.)*usvw2d[j][k]+6.*usvw2d[j][k+1]-3.*usvw2d[j][k+2]-(2./3.)*usvw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(wsuw2d[j][k+3]-wsuw2d[j][k+1])-(1./12.)*(wsuw2d[j][k+4]-wsuw2d[j][k]);
				Tuws2d[j][k] += 2.*nu*(-(7./3.)*wsuw2d[j][k]+6.*wsuw2d[j][k+1]-3.*wsuw2d[j][k+2]-(2./3.)*wsuw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(usww2d[j][k+3]-usww2d[j][k+1])-(1./12.)*(usww2d[j][k+4]-usww2d[j][k]);
				Tuws2d[j][k] += 2.*nu*(-(7./3.)*usww2d[j][k]+6.*usww2d[j][k+1]-3.*usww2d[j][k+2]-(2./3.)*usww2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(wswv2d[j][k+3]-wswv2d[j][k+1])-(1./12.)*(wswv2d[j][k+4]-wswv2d[j][k]);
				Tvws2d[j][k] += 2.*nu*(-(7./3.)*wswv2d[j][k]+6.*wswv2d[j][k+1]-3.*wswv2d[j][k+2]-(2./3.)*wswv2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(vsww2d[j][k+3]-vsww2d[j][k+1])-(1./12.)*(vsww2d[j][k+4]-vsww2d[j][k]);
				Tvws2d[j][k] += 2.*nu*(-(7./3.)*vsww2d[j][k]+6.*vsww2d[j][k+1]-3.*vsww2d[j][k+2]-(2./3.)*vsww2d[j][k+3]+3.*tmp);

				tmp = (2./3.)*(usuw2d[j][k+3]-usuw2d[j][k+1])-(1./12.)*(usuw2d[j][k+4]-usuw2d[j][k]);
				dudzusuw[j][k] = 2.*nu*(-(7./3.)*usuw2d[j][k]+6.*usuw2d[j][k+1]-3.*usuw2d[j][k+2]-(2./3.)*usuw2d[j][k+3]+3.*tmp);
				
				tmp = (2./3.)*(wke2d[j][k+3]-wke2d[j][k+1])-(1./12.)*(wke2d[j][k+4]-wke2d[j][k]);
				Tiir2d[j][k] -= ((-7./3.)*wke2d[j][k]+6.*wke2d[j][k+1]-3.*wke2d[j][k+2]-(2./3.)*wke2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(uuw2d[j][k+3]-uuw2d[j][k+1])-(1./12.)*(uuw2d[j][k+4]-uuw2d[j][k]);
				Tuur2d[j][k] -= ((-7./3)*uuw2d[j][k]+6.*uuw2d[j][k+1]-3.*uuw2d[j][k+2]-(2./3.)*uuw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(vvw2d[j][k+3]-vvw2d[j][k+1])-(1./12.)*(vvw2d[j][k+4]-vvw2d[j][k]);
				Tvvr2d[j][k] -= ((-7./3)*vvw2d[j][k]+6.*vvw2d[j][k+1]-3.*vvw2d[j][k+2]-(2./3.)*vvw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(www2d[j][k+3]-www2d[j][k+1])-(1./12.)*(www2d[j][k+4]-www2d[j][k]);
				Twwr2d[j][k] -= ((-7./3)*www2d[j][k]+6.*www2d[j][k+1]-3.*www2d[j][k+2]-(2./3.)*www2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(uvw2d[j][k+3]-uvw2d[j][k+1])-(1./12.)*(uvw2d[j][k+4]-uvw2d[j][k]);
				Tuvr2d[j][k] -= ((-7./3)*uvw2d[j][k]+6.*uvw2d[j][k+1]-3.*uvw2d[j][k+2]-(2./3.)*uvw2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(uww2d[j][k+3]-uww2d[j][k+1])-(1./12.)*(uww2d[j][k+4]-uww2d[j][k]);
				Tuwr2d[j][k] -= ((-7./3)*uww2d[j][k]+6.*uww2d[j][k+1]-3.*uww2d[j][k+2]-(2./3.)*uww2d[j][k+3]+3.*tmp);
				tmp = (2./3.)*(vww2d[j][k+3]-vww2d[j][k+1])-(1./12.)*(vww2d[j][k+4]-vww2d[j][k]);
				Tvwr2d[j][k] -= ((-7./3)*vww2d[j][k]+6.*vww2d[j][k+1]-3.*vww2d[j][k+2]-(2./3.)*vww2d[j][k+3]+3.*tmp);
			}
			else if (k==(zd-2))
			{
				tmp = (2./3.)*(uu2d[j][k-1]-uu2d[j][k-3])-(1./12.)*(uu2d[j][k]-uu2d[j][k-4]);
				Wduudz2d[j][k] = wAVG[j][k]*((7./3.)*uu2d[j][k]-6.*uu2d[j][k-1]+3.*uu2d[j][k-2]+(2./3.)*uu2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(vv2d[j][k-1]-vv2d[j][k-3])-(1./12.)*(vv2d[j][k]-vv2d[j][k-4]);
				Wdvvdz2d[j][k] = wAVG[j][k]*((7./3.)*vv2d[j][k]-6.*vv2d[j][k-1]+3.*vv2d[j][k-2]+(2./3.)*vv2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(ww2d[j][k-1]-ww2d[j][k-3])-(1./12.)*(ww2d[j][k]-ww2d[j][k-4]);
				Wdwwdz2d[j][k] = wAVG[j][k]*((7./3.)*ww2d[j][k]-6.*ww2d[j][k-1]+3.*ww2d[j][k-2]+(2./3.)*ww2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(uv2d[j][k-1]-uv2d[j][k-3])-(1./12.)*(uv2d[j][k]-uv2d[j][k-4]);
				Wduvdz2d[j][k] = wAVG[j][k]*((7./3.)*uv2d[j][k]-6.*uv2d[j][k-1]+3.*uv2d[j][k-2]+(2./3.)*uv2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(uw2d[j][k-1]-uw2d[j][k-3])-(1./12.)*(uw2d[j][k]-uw2d[j][k-4]);
				Wduwdz2d[j][k] = wAVG[j][k]*((7./3.)*uw2d[j][k]-6.*uw2d[j][k-1]+3.*uw2d[j][k-2]+(2./3.)*uw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(vw2d[j][k-1]-vw2d[j][k-3])-(1./12.)*(vw2d[j][k]-vw2d[j][k-4]);
				Wdvwdz2d[j][k] = wAVG[j][k]*((7./3.)*vw2d[j][k]-6.*vw2d[j][k-1]+3.*vw2d[j][k-2]+(2./3.)*vw2d[j][k-3]+3.*tmp);

				tmp = (2./3.)*(KE2d[j][k-1]-KE2d[j][k-3])-(1./12.)*(KE2d[j][k]-KE2d[j][k-4]);
				WdKdz2d[j][k] = wAVG[j][k]*((7./3.)*KE2d[j][k]-6.*KE2d[j][k-1]+3.*KE2d[j][k-2]+(2./3.)*KE2d[j][k-3]+3.*tmp);
				
				tmp = (2./3.)*(pw2d[j][k-1]-pw2d[j][k-3])-(1./12.)*(pw2d[j][k]-pw2d[j][k-4]);
				Tiipai2d[j][k] -= 0.5*2.*((7./3.)*pw2d[j][k]-6.*pw2d[j][k-1]+3.*pw2d[j][k-2]+(2./3.)*pw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(pw2d[j][k-1]-pw2d[j][k-3])-(1./12.)*(pw2d[j][k]-pw2d[j][k-4]);
				Twwpai2d[j][k] = -2.*((7./3.)*pw2d[j][k]-6.*pw2d[j][k-1]+3.*pw2d[j][k-2]+(2./3.)*pw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(pu2d[j][k-1]-pu2d[j][k-3])-(1./12.)*(pu2d[j][k]-pu2d[j][k-4]);
				Tuwpai2d[j][k] = -((7./3.)*pu2d[j][k]-6.*pu2d[j][k-1]+3.*pu2d[j][k-2]+(2./3.)*pu2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(pv2d[j][k-1]-pv2d[j][k-3])-(1./12.)*(pv2d[j][k]-pv2d[j][k-4]);
				Tvwpai2d[j][k] -= ((7./3.)*pv2d[j][k]-6.*pv2d[j][k-1]+3.*pv2d[j][k-2]+(2./3.)*pv2d[j][k-3]+3.*tmp);
				
				tmp = (2./3.)*(usuw2d[j][k-1]-usuw2d[j][k-3])-(1./12.)*(usuw2d[j][k]-usuw2d[j][k-4]);
				Tiis2d[j][k] += 2.*nu*((7./3.)*usuw2d[j][k]-6.*usuw2d[j][k-1]+3.*usuw2d[j][k-2]+(2./3.)*usuw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(vsvw2d[j][k-1]-vsvw2d[j][k-3])-(1./12.)*(vsvw2d[j][k]-vsvw2d[j][k-4]);
				Tiis2d[j][k] += 2.*nu*((7./3.)*vsvw2d[j][k]-6.*vsvw2d[j][k-1]+3.*vsvw2d[j][k-2]+(2./3.)*vsvw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(wsww2d[j][k-1]-wsww2d[j][k-3])-(1./12.)*(wsww2d[j][k]-wsww2d[j][k-4]);
				Tiis2d[j][k] += 2.*nu*((7./3.)*wsww2d[j][k]-6.*wsww2d[j][k-1]+3.*wsww2d[j][k-2]+(2./3.)*wsww2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(usuw2d[j][k-1]-usuw2d[j][k-3])-(1./12.)*(usuw2d[j][k]-usuw2d[j][k-4]);
				Tuus2d[j][k] += 4.*nu*((7./3.)*usuw2d[j][k]-6.*usuw2d[j][k-1]+3.*usuw2d[j][k-2]+(2./3.)*usuw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(vsvw2d[j][k-1]-vsvw2d[j][k-3])-(1./12.)*(vsvw2d[j][k]-vsvw2d[j][k-4]);
				Tvvs2d[j][k] += 4.*nu*((7./3.)*vsvw2d[j][k]-6.*vsvw2d[j][k-1]+3.*vsvw2d[j][k-2]+(2./3.)*vsvw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(wsww2d[j][k-1]-wsww2d[j][k-3])-(1./12.)*(wsww2d[j][k]-wsww2d[j][k-4]);
				Twws2d[j][k] += 4.*nu*((7./3.)*wsww2d[j][k]-6.*wsww2d[j][k-1]+3.*wsww2d[j][k-2]+(2./3.)*wsww2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(vsuw2d[j][k-1]-vsuw2d[j][k-3])-(1./12.)*(vsuw2d[j][k]-vsuw2d[j][k-4]);
				Tuvs2d[j][k] += 2.*nu*((7./3.)*vsuw2d[j][k]-6.*vsuw2d[j][k-1]+3.*vsuw2d[j][k-2]+(2./3.)*vsuw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(usvw2d[j][k-1]-usvw2d[j][k-3])-(1./12.)*(usvw2d[j][k]-usvw2d[j][k-4]);
				Tuvs2d[j][k] += 2.*nu*((7./3.)*usvw2d[j][k]-6.*usvw2d[j][k-1]+3.*usvw2d[j][k-2]+(2./3.)*usvw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(wsuw2d[j][k-1]-wsuw2d[j][k-3])-(1./12.)*(wsuw2d[j][k]-wsuw2d[j][k-4]);
				Tuws2d[j][k] += 2.*nu*((7./3.)*wsuw2d[j][k]-6.*wsuw2d[j][k-1]+3.*wsuw2d[j][k-2]+(2./3.)*wsuw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(usww2d[j][k-1]-usww2d[j][k-3])-(1./12.)*(usww2d[j][k]-usww2d[j][k-4]);
				Tuws2d[j][k] += 2.*nu*((7./3.)*usww2d[j][k]-6.*usww2d[j][k-1]+3.*usww2d[j][k-2]+(2./3.)*usww2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(wswv2d[j][k-1]-wswv2d[j][k-3])-(1./12.)*(wswv2d[j][k]-wswv2d[j][k-4]);
				Tvws2d[j][k] += 2.*nu*((7./3.)*wswv2d[j][k]-6.*wswv2d[j][k-1]+3.*wswv2d[j][k-2]+(2./3.)*wswv2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(vsww2d[j][k-1]-vsww2d[j][k-3])-(1./12.)*(vsww2d[j][k]-vsww2d[j][k-4]);
				Tvws2d[j][k] += 2.*nu*((7./3.)*vsww2d[j][k]-6.*vsww2d[j][k-1]+3.*vsww2d[j][k-2]+(2./3.)*vsww2d[j][k-3]+3.*tmp);

				tmp = (2./3.)*(usuw2d[j][k-1]-usuw2d[j][k-3])-(1./12.)*(usuw2d[j][k]-usuw2d[j][k-4]);
                                dudzusuw[j][k] = 2.*nu*((7./3.)*usuw2d[j][k]-6.*usuw2d[j][k-1]+3.*usuw2d[j][k-2]+(2./3.)*usuw2d[j][k-3]+3.*tmp);
				
				tmp = (2./3.)*(wke2d[j][k-1]-wke2d[j][k-3])-(1./12.)*(wke2d[j][k]-wke2d[j][k-4]);
				Tiir2d[j][k] -= ((7./3.)*wke2d[j][k]-6.*wke2d[j][k-1]+3.*wke2d[j][k-2]+(2./3.)*wke2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(uuw2d[j][k-1]-uuw2d[j][k-3])-(1./12.)*(uuw2d[j][k]-uuw2d[j][k-4]);
				Tuur2d[j][k] -= ((7./3.)*uuw2d[j][k]-6.*uuw2d[j][k-1]+3.*uuw2d[j][k-2]+(2./3.)*uuw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(vvw2d[j][k-1]-vvw2d[j][k-3])-(1./12.)*(vvw2d[j][k]-vvw2d[j][k-4]);
				Tvvr2d[j][k] -= ((7./3.)*vvw2d[j][k]-6.*vvw2d[j][k-1]+3.*vvw2d[j][k-2]+(2./3.)*vvw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(www2d[j][k-1]-www2d[j][k-3])-(1./12.)*(www2d[j][k]-www2d[j][k-4]);
				Twwr2d[j][k] -= ((7./3.)*www2d[j][k]-6.*www2d[j][k-1]+3.*www2d[j][k-2]+(2./3.)*www2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(uvw2d[j][k-1]-uvw2d[j][k-3])-(1./12.)*(uvw2d[j][k]-uvw2d[j][k-4]);
				Tuvr2d[j][k] -= ((7./3.)*uvw2d[j][k]-6.*uvw2d[j][k-1]+3.*uvw2d[j][k-2]+(2./3.)*uvw2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(uww2d[j][k-1]-uww2d[j][k-3])-(1./12.)*(uww2d[j][k]-uww2d[j][k-4]);
				Tuwr2d[j][k] -= ((7./3.)*uww2d[j][k]-6.*uww2d[j][k-1]+3.*uww2d[j][k-2]+(2./3.)*uww2d[j][k-3]+3.*tmp);
				tmp = (2./3.)*(vww2d[j][k-1]-vww2d[j][k-3])-(1./12.)*(vww2d[j][k]-vww2d[j][k-4]);
				Tvwr2d[j][k] -= ((7./3.)*vww2d[j][k]-6.*vww2d[j][k-1]+3.*vww2d[j][k-2]+(2./3.)*vww2d[j][k-3]+3.*tmp);
			}
			else if (k==2)
			{
				tmp = (2./3.)*(uu2d[j][k+2]-uu2d[j][k])-(1./12.)*(uu2d[j][k+3]-uu2d[j][k-1]);
				Wduudz2d[j][k] = wAVG[j][k]*(-(1./6)*uu2d[j][k-1]-(3./2.)*uu2d[j][k]+(3./2.)*uu2d[j][k+1]+(1./6.)*uu2d[j][k+2]-tmp);
				tmp = (2./3.)*(vv2d[j][k+2]-vv2d[j][k])-(1./12.)*(vv2d[j][k+3]-vv2d[j][k-1]);
				Wdvvdz2d[j][k] = wAVG[j][k]*(-(1./6)*vv2d[j][k-1]-(3./2.)*vv2d[j][k]+(3./2.)*vv2d[j][k+1]+(1./6.)*vv2d[j][k+2]-tmp);
				tmp = (2./3.)*(ww2d[j][k+2]-ww2d[j][k])-(1./12.)*(ww2d[j][k+3]-ww2d[j][k-1]);
				Wdwwdz2d[j][k] = wAVG[j][k]*(-(1./6)*ww2d[j][k-1]-(3./2.)*ww2d[j][k]+(3./2.)*ww2d[j][k+1]+(1./6.)*ww2d[j][k+2]-tmp);
				tmp = (2./3.)*(uv2d[j][k+2]-uv2d[j][k])-(1./12.)*(uv2d[j][k+3]-uv2d[j][k-1]);
				Wduvdz2d[j][k] = wAVG[j][k]*(-(1./6)*uv2d[j][k-1]-(3./2.)*uv2d[j][k]+(3./2.)*uv2d[j][k+1]+(1./6.)*uv2d[j][k+2]-tmp);
				tmp = (2./3.)*(uw2d[j][k+2]-uw2d[j][k])-(1./12.)*(uw2d[j][k+3]-uw2d[j][k-1]);
				Wduwdz2d[j][k] = wAVG[j][k]*(-(1./6)*uw2d[j][k-1]-(3./2.)*uw2d[j][k]+(3./2.)*uw2d[j][k+1]+(1./6.)*uw2d[j][k+2]-tmp);
				tmp = (2./3.)*(vw2d[j][k+2]-vw2d[j][k])-(1./12.)*(vw2d[j][k+3]-vw2d[j][k-1]);
				Wdvwdz2d[j][k] = wAVG[j][k]*(-(1./6)*vw2d[j][k-1]-(3./2.)*vw2d[j][k]+(3./2.)*vw2d[j][k+1]+(1./6.)*vw2d[j][k+2]-tmp);

				tmp = (2./3.)*(KE2d[j][k+2]-KE2d[j][k])-(1./12.)*(KE2d[j][k+3]-KE2d[j][k-1]);
				WdKdz2d[j][k] = wAVG[j][k]*(-(1./6)*KE2d[j][k-1]-(3./2.)*KE2d[j][k]+(3./2.)*KE2d[j][k+1]+(1./6.)*KE2d[j][k+2]-tmp);
				
				tmp = (2./3.)*(pw2d[j][k+2]-pw2d[j][k])-(1./12.)*(pw2d[j][k+3]-pw2d[j][k-1]);
				Tiipai2d[j][k] -= 0.5*2.*(-(1./6)*pw2d[j][k-1]-(3./2.)*pw2d[j][k]+(3./2.)*pw2d[j][k+1]+(1./6.)*pw2d[j][k+2]-tmp);
				tmp = (2./3.)*(pw2d[j][k+2]-pw2d[j][k])-(1./12.)*(pw2d[j][k+3]-pw2d[j][k-1]);
				Twwpai2d[j][k] = -2.*(-(1./6)*pw2d[j][k-1]-(3./2.)*pw2d[j][k]+(3./2.)*pw2d[j][k+1]+(1./6.)*pw2d[j][k+2]-tmp);
				tmp = (2./3.)*(pu2d[j][k+2]-pu2d[j][k])-(1./12.)*(pu2d[j][k+3]-pu2d[j][k-1]);
				Tuwpai2d[j][k] = (-(1./6)*pu2d[j][k-1]-(3./2.)*pu2d[j][k]+(3./2.)*pu2d[j][k+1]+(1./6.)*pu2d[j][k+2]-tmp);
				tmp = (2./3.)*(pv2d[j][k+2]-pv2d[j][k])-(1./12.)*(pv2d[j][k+3]-pv2d[j][k-1]);
				Tvwpai2d[j][k] -= (-(1./6)*pv2d[j][k-1]-(3./2.)*pv2d[j][k]+(3./2.)*pv2d[j][k+1]+(1./6.)*pv2d[j][k+2]-tmp);
				
				tmp = (2./3.)*(usuw2d[j][k+2]-usuw2d[j][k])-(1./12.)*(usuw2d[j][k+3]-usuw2d[j][k-1]);
				Tiis2d[j][k] += 2.*nu*(-(1./6)*usuw2d[j][k-1]-(3./2.)*usuw2d[j][k]+(3./2.)*usuw2d[j][k+1]+(1./6.)*usuw2d[j][k+2]-tmp);
				tmp = (2./3.)*(vsvw2d[j][k+2]-vsvw2d[j][k])-(1./12.)*(vsvw2d[j][k+3]-vsvw2d[j][k-1]);
				Tiis2d[j][k] += 2.*nu*(-(1./6)*vsvw2d[j][k-1]-(3./2.)*vsvw2d[j][k]+(3./2.)*vsvw2d[j][k+1]+(1./6.)*vsvw2d[j][k+2]-tmp);
				tmp = (2./3.)*(wsww2d[j][k+2]-wsww2d[j][k])-(1./12.)*(wsww2d[j][k+3]-wsww2d[j][k-1]);
				Tiis2d[j][k] += 2.*nu*(-(1./6)*wsww2d[j][k-1]-(3./2.)*wsww2d[j][k]+(3./2.)*wsww2d[j][k+1]+(1./6.)*wsww2d[j][k+2]-tmp);
				tmp = (2./3.)*(usuw2d[j][k+2]-usuw2d[j][k])-(1./12.)*(usuw2d[j][k+3]-usuw2d[j][k-1]);
				Tuus2d[j][k] += 4.*nu*(-(1./6)*usuw2d[j][k-1]-(3./2.)*usuw2d[j][k]+(3./2.)*usuw2d[j][k+1]+(1./6.)*usuw2d[j][k+2]-tmp);
				tmp = (2./3.)*(vsvw2d[j][k+2]-vsvw2d[j][k])-(1./12.)*(vsvw2d[j][k+3]-vsvw2d[j][k-1]);
				Tvvs2d[j][k] += 4.*nu*(-(1./6)*vsvw2d[j][k-1]-(3./2.)*vsvw2d[j][k]+(3./2.)*vsvw2d[j][k+1]+(1./6.)*vsvw2d[j][k+2]-tmp);
				tmp = (2./3.)*(wsww2d[j][k+2]-wsww2d[j][k])-(1./12.)*(wsww2d[j][k+3]-wsww2d[j][k-1]);
				Twws2d[j][k] += 4.*nu*(-(1./6)*wsww2d[j][k-1]-(3./2.)*wsww2d[j][k]+(3./2.)*wsww2d[j][k+1]+(1./6.)*wsww2d[j][k+2]-tmp);
				tmp = (2./3.)*(vsuw2d[j][k+2]-vsuw2d[j][k])-(1./12.)*(vsuw2d[j][k+3]-vsuw2d[j][k-1]);
				Tuvs2d[j][k] += 2.*nu*(-(1./6)*vsuw2d[j][k-1]-(3./2.)*vsuw2d[j][k]+(3./2.)*vsuw2d[j][k+1]+(1./6.)*vsuw2d[j][k+2]-tmp);
				tmp = (2./3.)*(usvw2d[j][k+2]-usvw2d[j][k])-(1./12.)*(usvw2d[j][k+3]-usvw2d[j][k-1]);
				Tuvs2d[j][k] += 2.*nu*(-(1./6)*usvw2d[j][k-1]-(3./2.)*usvw2d[j][k]+(3./2.)*usvw2d[j][k+1]+(1./6.)*usvw2d[j][k+2]-tmp);
				tmp = (2./3.)*(wsuw2d[j][k+2]-wsuw2d[j][k])-(1./12.)*(wsuw2d[j][k+3]-wsuw2d[j][k-1]);
				Tuws2d[j][k] += 2.*nu*(-(1./6)*wsuw2d[j][k-1]-(3./2.)*wsuw2d[j][k]+(3./2.)*wsuw2d[j][k+1]+(1./6.)*wsuw2d[j][k+2]-tmp);
				tmp = (2./3.)*(usww2d[j][k+2]-usww2d[j][k])-(1./12.)*(usww2d[j][k+3]-usww2d[j][k-1]);
				Tuws2d[j][k] += 2.*nu*(-(1./6)*usww2d[j][k-1]-(3./2.)*usww2d[j][k]+(3./2.)*usww2d[j][k+1]+(1./6.)*usww2d[j][k+2]-tmp);
				tmp = (2./3.)*(wswv2d[j][k+2]-wswv2d[j][k])-(1./12.)*(wswv2d[j][k+3]-wswv2d[j][k-1]);
				Tvws2d[j][k] += 2.*nu*(-(1./6)*wswv2d[j][k-1]-(3./2.)*wswv2d[j][k]+(3./2.)*wswv2d[j][k+1]+(1./6.)*wswv2d[j][k+2]-tmp);
				tmp = (2./3.)*(vsww2d[j][k+2]-vsww2d[j][k])-(1./12.)*(vsww2d[j][k+3]-vsww2d[j][k-1]);
				Tvws2d[j][k] += 2.*nu*(-(1./6)*vsww2d[j][k-1]-(3./2.)*vsww2d[j][k]+(3./2.)*vsww2d[j][k+1]+(1./6.)*vsww2d[j][k+2]-tmp);
				
				tmp = (2./3.)*(usuw2d[j][k+2]-usuw2d[j][k])-(1./12.)*(usuw2d[j][k+3]-usuw2d[j][k-1]);
                                dudzusuw[j][k] = 2.*nu*(-(1./6)*usuw2d[j][k-1]-(3./2.)*usuw2d[j][k]+(3./2.)*usuw2d[j][k+1]+(1./6.)*usuw2d[j][k+2]-tmp);

				tmp = (2./3.)*(wke2d[j][k+2]-wke2d[j][k])-(1./12.)*(wke2d[j][k+3]-wke2d[j][k-1]);
				Tiir2d[j][k] -= (-(1./6)*wke2d[j][k-1]-(3./2.)*wke2d[j][k]+(3./2.)*wke2d[j][k+1]+(1./6.)*wke2d[j][k+2]-tmp);
				tmp = (2./3.)*(uuw2d[j][k+2]-uuw2d[j][k])-(1./12.)*(uuw2d[j][k+3]-uuw2d[j][k-1]);
				Tuur2d[j][k] -= (-(1./6)*uuw2d[j][k-1]-(3./2.)*uuw2d[j][k]+(3./2.)*uuw2d[j][k+1]+(1./6.)*uuw2d[j][k+2]-tmp);
				tmp = (2./3.)*(vvw2d[j][k+2]-vvw2d[j][k])-(1./12.)*(vvw2d[j][k+3]-vvw2d[j][k-1]);
				Tvvr2d[j][k] -= (-(1./6)*vvw2d[j][k-1]-(3./2.)*vvw2d[j][k]+(3./2.)*vvw2d[j][k+1]+(1./6.)*vvw2d[j][k+2]-tmp);
				tmp = (2./3.)*(www2d[j][k+2]-www2d[j][k])-(1./12.)*(www2d[j][k+3]-www2d[j][k-1]);
				Twwr2d[j][k] -= (-(1./6)*www2d[j][k-1]-(3./2.)*www2d[j][k]+(3./2.)*www2d[j][k+1]+(1./6.)*www2d[j][k+2]-tmp);
				tmp = (2./3.)*(uvw2d[j][k+2]-uvw2d[j][k])-(1./12.)*(uvw2d[j][k+3]-uvw2d[j][k-1]);
				Tuvr2d[j][k] -= (-(1./6)*uvw2d[j][k-1]-(3./2.)*uvw2d[j][k]+(3./2.)*uvw2d[j][k+1]+(1./6.)*uvw2d[j][k+2]-tmp);
				tmp = (2./3.)*(uww2d[j][k+2]-uww2d[j][k])-(1./12.)*(uww2d[j][k+3]-uww2d[j][k-1]);
				Tuwr2d[j][k] -= (-(1./6)*uww2d[j][k-1]-(3./2.)*uww2d[j][k]+(3./2.)*uww2d[j][k+1]+(1./6.)*uww2d[j][k+2]-tmp);
				tmp = (2./3.)*(vww2d[j][k+2]-vww2d[j][k])-(1./12.)*(vww2d[j][k+3]-vww2d[j][k-1]);
				Tvwr2d[j][k] -= (-(1./6)*vww2d[j][k-1]-(3./2.)*vww2d[j][k]+(3./2.)*vww2d[j][k+1]+(1./6.)*vww2d[j][k+2]-tmp);
			}
			else if (k==(zd-3))
			{
				tmp = (2./3.)*(uu2d[j][k]-uu2d[j][k-2])-(1./12.)*(uu2d[j][k+1]-uu2d[j][k-3]);
				Wduudz2d[j][k] = wAVG[j][k]*((1./6.)*uu2d[j][k+1]+(3./2.)*uu2d[j][k]-(3./2.)*uu2d[j][k-1]-(1./6.)*uu2d[j][k-2]-tmp);
				tmp = (2./3.)*(vv2d[j][k]-vv2d[j][k-2])-(1./12.)*(vv2d[j][k+1]-vv2d[j][k-3]);
				Wdvvdz2d[j][k] = wAVG[j][k]*((1./6.)*vv2d[j][k+1]+(3./2.)*vv2d[j][k]-(3./2.)*vv2d[j][k-1]-(1./6.)*vv2d[j][k-2]-tmp);
				tmp = (2./3.)*(ww2d[j][k]-ww2d[j][k-2])-(1./12.)*(ww2d[j][k+1]-ww2d[j][k-3]);
				Wdwwdz2d[j][k] = wAVG[j][k]*((1./6.)*ww2d[j][k+1]+(3./2.)*ww2d[j][k]-(3./2.)*ww2d[j][k-1]-(1./6.)*ww2d[j][k-2]-tmp);
				tmp = (2./3.)*(uv2d[j][k]-uv2d[j][k-2])-(1./12.)*(uv2d[j][k+1]-uv2d[j][k-3]);
				Wduvdz2d[j][k] = wAVG[j][k]*((1./6.)*uv2d[j][k+1]+(3./2.)*uv2d[j][k]-(3./2.)*uv2d[j][k-1]-(1./6.)*uv2d[j][k-2]-tmp);
				tmp = (2./3.)*(uw2d[j][k]-uw2d[j][k-2])-(1./12.)*(uw2d[j][k+1]-uw2d[j][k-3]);
				Wduwdz2d[j][k] = wAVG[j][k]*((1./6.)*uw2d[j][k+1]+(3./2.)*uw2d[j][k]-(3./2.)*uw2d[j][k-1]-(1./6.)*uw2d[j][k-2]-tmp);
				tmp = (2./3.)*(vw2d[j][k]-vw2d[j][k-2])-(1./12.)*(vw2d[j][k+1]-vw2d[j][k-3]);
				Wdvwdz2d[j][k] = wAVG[j][k]*((1./6.)*vw2d[j][k+1]+(3./2.)*vw2d[j][k]-(3./2.)*vw2d[j][k-1]-(1./6.)*vw2d[j][k-2]-tmp);

				tmp = (2./3.)*(KE2d[j][k]-KE2d[j][k-2])-(1./12.)*(KE2d[j][k+1]-KE2d[j][k-3]);
				WdKdz2d[j][k] = wAVG[j][k]*((1./6.)*KE2d[j][k+1]+(3./2.)*KE2d[j][k]-(3./2.)*KE2d[j][k-1]-(1./6.)*KE2d[j][k-2]-tmp);
				
				tmp = (2./3.)*(pw2d[j][k]-pw2d[j][k-2])-(1./12.)*(pw2d[j][k+1]-pw2d[j][k-3]);
				Tiipai2d[j][k] -= 0.5*2.*((1./6.)*pw2d[j][k+1]+(3./2.)*pw2d[j][k]-(3./2.)*pw2d[j][k-1]-(1./6.)*pw2d[j][k-2]-tmp);
				tmp = (2./3.)*(pw2d[j][k]-pw2d[j][k-2])-(1./12.)*(pw2d[j][k+1]-pw2d[j][k-3]);
				Twwpai2d[j][k] = -2.*((1./6.)*pw2d[j][k+1]+(3./2.)*pw2d[j][k]-(3./2.)*pw2d[j][k-1]-(1./6.)*pw2d[j][k-2]-tmp);
				tmp = (2./3.)*(pu2d[j][k]-pu2d[j][k-2])-(1./12.)*(pu2d[j][k+1]-pu2d[j][k-3]);
				Tuwpai2d[j][k] = -((1./6.)*pu2d[j][k+1]+(3./2.)*pu2d[j][k]-(3./2.)*pu2d[j][k-1]-(1./6.)*pu2d[j][k-2]-tmp);
				tmp = (2./3.)*(pv2d[j][k]-pv2d[j][k-2])-(1./12.)*(pv2d[j][k+1]-pv2d[j][k-3]);
				Tvwpai2d[j][k] -= ((1./6.)*pv2d[j][k+1]+(3./2.)*pv2d[j][k]-(3./2.)*pv2d[j][k-1]-(1./6.)*pv2d[j][k-2]-tmp);
				
				tmp = (2./3.)*(usuw2d[j][k]-usuw2d[j][k-2])-(1./12.)*(usuw2d[j][k+1]-usuw2d[j][k-3]);
				Tiis2d[j][k] += 2.*nu*((1./6.)*usuw2d[j][k+1]+(3./2.)*usuw2d[j][k]-(3./2.)*usuw2d[j][k-1]-(1./6.)*usuw2d[j][k-2]-tmp);
				tmp = (2./3.)*(vsvw2d[j][k]-vsvw2d[j][k-2])-(1./12.)*(vsvw2d[j][k+1]-vsvw2d[j][k-3]);
				Tiis2d[j][k] += 2.*nu*((1./6.)*vsvw2d[j][k+1]+(3./2.)*vsvw2d[j][k]-(3./2.)*vsvw2d[j][k-1]-(1./6.)*vsvw2d[j][k-2]-tmp);
				tmp = (2./3.)*(wsww2d[j][k]-wsww2d[j][k-2])-(1./12.)*(wsww2d[j][k+1]-wsww2d[j][k-3]);
				Tiis2d[j][k] += 2.*nu*((1./6.)*wsww2d[j][k+1]+(3./2.)*wsww2d[j][k]-(3./2.)*wsww2d[j][k-1]-(1./6.)*wsww2d[j][k-2]-tmp);
				tmp = (2./3.)*(usuw2d[j][k]-usuw2d[j][k-2])-(1./12.)*(usuw2d[j][k+1]-usuw2d[j][k-3]);
				Tuus2d[j][k] += 4.*nu*((1./6.)*usuw2d[j][k+1]+(3./2.)*usuw2d[j][k]-(3./2.)*usuw2d[j][k-1]-(1./6.)*usuw2d[j][k-2]-tmp);
				tmp = (2./3.)*(vsvw2d[j][k]-vsvw2d[j][k-2])-(1./12.)*(vsvw2d[j][k+1]-vsvw2d[j][k-3]);
				Tvvs2d[j][k] += 4.*nu*((1./6.)*vsvw2d[j][k+1]+(3./2.)*vsvw2d[j][k]-(3./2.)*vsvw2d[j][k-1]-(1./6.)*vsvw2d[j][k-2]-tmp);
				tmp = (2./3.)*(wsww2d[j][k]-wsww2d[j][k-2])-(1./12.)*(wsww2d[j][k+1]-wsww2d[j][k-3]);
				Twws2d[j][k] += 4.*nu*((1./6.)*wsww2d[j][k+1]+(3./2.)*wsww2d[j][k]-(3./2.)*wsww2d[j][k-1]-(1./6.)*wsww2d[j][k-2]-tmp);
				tmp = (2./3.)*(vsuw2d[j][k]-vsuw2d[j][k-2])-(1./12.)*(vsuw2d[j][k+1]-vsuw2d[j][k-3]);
				Tuvs2d[j][k] += 2.*nu*((1./6.)*vsuw2d[j][k+1]+(3./2.)*vsuw2d[j][k]-(3./2.)*vsuw2d[j][k-1]-(1./6.)*vsuw2d[j][k-2]-tmp);
				tmp = (2./3.)*(usvw2d[j][k]-usvw2d[j][k-2])-(1./12.)*(usvw2d[j][k+1]-usvw2d[j][k-3]);
				Tuvs2d[j][k] += 2.*nu*((1./6.)*usvw2d[j][k+1]+(3./2.)*usvw2d[j][k]-(3./2.)*usvw2d[j][k-1]-(1./6.)*usvw2d[j][k-2]-tmp);
				tmp = (2./3.)*(wsuw2d[j][k]-wsuw2d[j][k-2])-(1./12.)*(wsuw2d[j][k+1]-wsuw2d[j][k-3]);
				Tuws2d[j][k] += 2.*nu*((1./6.)*wsuw2d[j][k+1]+(3./2.)*wsuw2d[j][k]-(3./2.)*wsuw2d[j][k-1]-(1./6.)*wsuw2d[j][k-2]-tmp);
				tmp = (2./3.)*(usww2d[j][k]-usww2d[j][k-2])-(1./12.)*(usww2d[j][k+1]-usww2d[j][k-3]);
				Tuws2d[j][k] += 2.*nu*((1./6.)*usww2d[j][k+1]+(3./2.)*usww2d[j][k]-(3./2.)*usww2d[j][k-1]-(1./6.)*usww2d[j][k-2]-tmp);
				tmp = (2./3.)*(wswv2d[j][k]-wswv2d[j][k-2])-(1./12.)*(wswv2d[j][k+1]-wswv2d[j][k-3]);
				Tvws2d[j][k] += 2.*nu*((1./6.)*wswv2d[j][k+1]+(3./2.)*wswv2d[j][k]-(3./2.)*wswv2d[j][k-1]-(1./6.)*wswv2d[j][k-2]-tmp);
				tmp = (2./3.)*(vsww2d[j][k]-vsww2d[j][k-2])-(1./12.)*(vsww2d[j][k+1]-vsww2d[j][k-3]);
				Tvws2d[j][k] += 2.*nu*((1./6.)*vsww2d[j][k+1]+(3./2.)*vsww2d[j][k]-(3./2.)*vsww2d[j][k-1]-(1./6.)*vsww2d[j][k-2]-tmp);
				
				tmp = (2./3.)*(usuw2d[j][k]-usuw2d[j][k-2])-(1./12.)*(usuw2d[j][k+1]-usuw2d[j][k-3]);
                                dudzusuw[j][k] = 2.*nu*((1./6.)*usuw2d[j][k+1]+(3./2.)*usuw2d[j][k]-(3./2.)*usuw2d[j][k-1]-(1./6.)*usuw2d[j][k-2]-tmp);

				tmp = (2./3.)*(wke2d[j][k]-wke2d[j][k-2])-(1./12.)*(wke2d[j][k+1]-wke2d[j][k-3]);
				Tiir2d[j][k] -= ((1./6.)*wke2d[j][k+1]+(3./2.)*wke2d[j][k]-(3./2.)*wke2d[j][k-1]-(1./6.)*wke2d[j][k-2]-tmp);
				tmp = (2./3.)*(uuw2d[j][k]-uuw2d[j][k-2])-(1./12.)*(uuw2d[j][k+1]-uuw2d[j][k-3]);
				Tuur2d[j][k] -= ((1./6.)*uuw2d[j][k+1]+(3./2.)*uuw2d[j][k]-(3./2.)*uuw2d[j][k-1]-(1./6.)*uuw2d[j][k-2]-tmp);
				tmp = (2./3.)*(vvw2d[j][k]-vvw2d[j][k-2])-(1./12.)*(vvw2d[j][k+1]-vvw2d[j][k-3]);
				Tvvr2d[j][k] -= ((1./6.)*vvw2d[j][k+1]+(3./2.)*vvw2d[j][k]-(3./2.)*vvw2d[j][k-1]-(1./6.)*vvw2d[j][k-2]-tmp);
				tmp = (2./3.)*(www2d[j][k]-www2d[j][k-2])-(1./12.)*(www2d[j][k+1]-www2d[j][k-3]);
				Twwr2d[j][k] -= ((1./6.)*www2d[j][k+1]+(3./2.)*www2d[j][k]-(3./2.)*www2d[j][k-1]-(1./6.)*www2d[j][k-2]-tmp);
				tmp = (2./3.)*(uvw2d[j][k]-uvw2d[j][k-2])-(1./12.)*(uvw2d[j][k+1]-uvw2d[j][k-3]);
				Tuvr2d[j][k] -= ((1./6.)*uvw2d[j][k+1]+(3./2.)*uvw2d[j][k]-(3./2.)*uvw2d[j][k-1]-(1./6.)*uvw2d[j][k-2]-tmp);
				tmp = (2./3.)*(uww2d[j][k]-uww2d[j][k-2])-(1./12.)*(uww2d[j][k+1]-uww2d[j][k-3]);
				Tuwr2d[j][k] -= ((1./6.)*uww2d[j][k+1]+(3./2.)*uww2d[j][k]-(3./2.)*uww2d[j][k-1]-(1./6.)*uww2d[j][k-2]-tmp);
				tmp = (2./3.)*(vww2d[j][k]-vww2d[j][k-2])-(1./12.)*(vww2d[j][k+1]-vww2d[j][k-3]);
				Tvwr2d[j][k] -= ((1./6.)*vww2d[j][k+1]+(3./2.)*vww2d[j][k]-(3./2.)*vww2d[j][k-1]-(1./6.)*vww2d[j][k-2]-tmp);
			}
			else
			{
				Wduudz2d[j][k] = wAVG[j][k]*((2./3.)*(uu2d[j][k+1]-uu2d[j][k-1])-(1./12.)*(uu2d[j][k+2]-uu2d[j][k-2]));
				Wdvvdz2d[j][k] = wAVG[j][k]*((2./3.)*(vv2d[j][k+1]-vv2d[j][k-1])-(1./12.)*(vv2d[j][k+2]-vv2d[j][k-2]));
				Wdwwdz2d[j][k] = wAVG[j][k]*((2./3.)*(ww2d[j][k+1]-ww2d[j][k-1])-(1./12.)*(ww2d[j][k+2]-ww2d[j][k-2]));
				Wduvdz2d[j][k] = wAVG[j][k]*((2./3.)*(uv2d[j][k+1]-uv2d[j][k-1])-(1./12.)*(uv2d[j][k+2]-uv2d[j][k-2]));
				Wduwdz2d[j][k] = wAVG[j][k]*((2./3.)*(uw2d[j][k+1]-uw2d[j][k-1])-(1./12.)*(uw2d[j][k+2]-uw2d[j][k-2]));
				Wdvwdz2d[j][k] = wAVG[j][k]*((2./3.)*(vw2d[j][k+1]-vw2d[j][k-1])-(1./12.)*(vw2d[j][k+2]-vw2d[j][k-2]));

				WdKdz2d[j][k] = wAVG[j][k]*((2./3.)*(KE2d[j][k+1]-KE2d[j][k-1])-(1./12.)*(KE2d[j][k+2]-KE2d[j][k-2]));
				
				Tiipai2d[j][k] -= 0.5*2.*((2./3.)*(pw2d[j][k+1]-pw2d[j][k-1])-(1./12.)*(pw2d[j][k+2]-pw2d[j][k-2]));
				Twwpai2d[j][k] = -2.*((2./3.)*(pw2d[j][k+1]-pw2d[j][k-1])-(1./12.)*(pw2d[j][k+2]-pw2d[j][k-2]));
				Tuwpai2d[j][k] = -((2./3.)*(pu2d[j][k+1]-pu2d[j][k-1])-(1./12.)*(pu2d[j][k+2]-pu2d[j][k-2]));
				Tvwpai2d[j][k] -= ((2./3.)*(pv2d[j][k+1]-pv2d[j][k-1])-(1./12.)*(pv2d[j][k+2]-pv2d[j][k-2]));
				
				Tiis2d[j][k] += 2.*nu*((2./3.)*(usuw2d[j][k+1]-usuw2d[j][k-1]+vsvw2d[j][k+1]-vsvw2d[j][k-1]+wsww2d[j][k+1]-wsww2d[j][k-1])-(1./12.)*(usuw2d[j][k+2]-usuw2d[j][k-2]+vsvw2d[j][k+2]-vsvw2d[j][k-2]+wsww2d[j][k+2]-wsww2d[j][k-2]));
				Tuus2d[j][k] += 4.*nu*((2./3.)*(usuw2d[j][k+1]-usuw2d[j][k-1])-(1./12.)*(usuw2d[j][k+2]-usuw2d[j][k-2]));
				Tvvs2d[j][k] += 4.*nu*((2./3.)*(vsvw2d[j][k+1]-vsvw2d[j][k-1])-(1./12.)*(vsvw2d[j][k+2]-vsvw2d[j][k-2]));
				Twws2d[j][k] += 4.*nu*((2./3.)*(wsww2d[j][k+1]-wsww2d[j][k-1])-(1./12.)*(wsww2d[j][k+2]-wsww2d[j][k-2]));
				Tuvs2d[j][k] += 2.*nu*((2./3.)*(vsuw2d[j][k+1]-vsuw2d[j][k-1]+usvw2d[j][k+1]-usvw2d[j][k-1])-(1./12.)*(vsuw2d[j][k+2]-vsuw2d[j][k-2]+usvw2d[j][k+2]-usvw2d[j][k-2]));
				Tuws2d[j][k] += 2.*nu*((2./3.)*(wsuw2d[j][k+1]-wsuw2d[j][k-1]+usww2d[j][k+1]-usww2d[j][k-1])-(1./12.)*(wsuw2d[j][k+2]-wsuw2d[j][k-2]+usww2d[j][k+2]-usww2d[j][k-2]));
				Tvws2d[j][k] += 2.*nu*((2./3.)*(wswv2d[j][k+1]-wswv2d[j][k-1]+vsww2d[j][k+1]-vsww2d[j][k-1])-(1./12.)*(wswv2d[j][k+2]-wswv2d[j][k-2]+vsww2d[j][k+2]-vsww2d[j][k-2]));
				
				dudzusuw[j][k] = 2.*nu*((2./3.)*(usuw2d[j][k+1]-usuw2d[j][k-1])-(1./12.)*(usuw2d[j][k+2]-usuw2d[j][k-2]));

				Tiir2d[j][k] -= ((2./3.)*(wke2d[j][k+1]-wke2d[j][k-1])-(1./12.)*(wke2d[j][k+2]-wke2d[j][k-2]));
				Tuur2d[j][k] -= ((2./3.)*(uuw2d[j][k+1]-uuw2d[j][k-1])-(1./12.)*(uuw2d[j][k+2]-uuw2d[j][k-2]));
				Tvvr2d[j][k] -= ((2./3.)*(vvw2d[j][k+1]-vvw2d[j][k-1])-(1./12.)*(vvw2d[j][k+2]-vvw2d[j][k-2]));
				Twwr2d[j][k] -= ((2./3.)*(www2d[j][k+1]-www2d[j][k-1])-(1./12.)*(www2d[j][k+2]-www2d[j][k-2]));
				Tuvr2d[j][k] -= ((2./3.)*(uvw2d[j][k+1]-uvw2d[j][k-1])-(1./12.)*(uvw2d[j][k+2]-uvw2d[j][k-2]));
				Tuwr2d[j][k] -= ((2./3.)*(uww2d[j][k+1]-uww2d[j][k-1])-(1./12.)*(uww2d[j][k+2]-uww2d[j][k-2]));
				Tvwr2d[j][k] -= ((2./3.)*(vww2d[j][k+1]-vww2d[j][k-1])-(1./12.)*(vww2d[j][k+2]-vww2d[j][k-2]));
			}
		}
	}
	
	for (k=0;k<zd;k++)
	{
		p1[k] = 0.;
		p2[k] = 0.;
		p3[k] = 0.;
		p4[k] = 0.;
		p5[k] = 0.;
		p6[k] = 0.;
		p7[k] = 0.;
		p8[k] = 0.;
		p9[k] = 0.;
		p10[k] = 0.;
		p11[k] = 0.;
		p12[k] = 0.;
		p13[k] = 0.;
		p14[k] = 0.;
		p15[k] = 0.;
		p16[k] = 0.;
		p17[k] = 0.;
		p18[k] = 0.;
	}

	for (k=1;k<(zd-1);k++)
	{
		siu = *(UTIN+k) = sqrt((*(uu+k)));
		siv = *(VTIN+k) = sqrt((*(vv+k)));
		siw = *(WTIN+k) = sqrt((*(ww+k)));
		sip = *(DTIN+k) = sqrt((*(dd+k)));
		
		*(wx+k) = sqrt((*(wx+k)));
		*(wy+k) = sqrt((*(wy+k)));
		*(wz+k) = sqrt((*(wz+k)));
		
		*(SKEWu+k) /= (siu*siu*siu);
		*(SKEWv+k) /= (siv*siv*siv);
		*(SKEWw+k) /= (siw*siw*siw);
		*(SKEWp+k) /= (sip*sip*sip);
		*(SKEWuw+k) /= (siu*siu*siu*siw*siw*siw);
		*(SKEWvw+k) /= (siv*siv*siv*siw*siw*siw);
		*(SKEWuv+k) /= (siu*siu*siu*siv*siv*siv);
		
		*(KURTu+k) /= (siu*siu*siu*siu);
		*(KURTv+k) /= (siv*siv*siv*siv);
		*(KURTw+k) /= (siw*siw*siw*siw);
		*(KURTp+k) /= (sip*sip*sip*sip);
		*(KURTuw+k) /= (siu*siu*siu*siu*siw*siw*siw*siw);
		*(KURTuv+k) /= (siu*siu*siu*siu*siv*siv*siv*siv);
		*(KURTvw+k) /= (siv*siv*siv*siv*siw*siw*siw*siw);
		
		for (j=0;j<exty;j++)
		{
			*(KE+k) += KE2d[j][k];
			
			*(Pii+k) += Pii2d[j][k];
			*(Puu+k) += Puu2d[j][k];
			*(Pvv+k) += Pvv2d[j][k];
			*(Pww+k) += Pww2d[j][k];
			*(Puv+k) += Puv2d[j][k];
			*(Puw+k) += Puw2d[j][k];
			*(Pvw+k) += Pvw2d[j][k];
			
			*(Tiipai+k) += Tiipai2d[j][k];
			*(Tuupai+k) += Tuupai2d[j][k];
			*(Tvvpai+k) += Tvvpai2d[j][k];
			*(Twwpai+k) += Twwpai2d[j][k];
			*(Tuvpai+k) += Tuvpai2d[j][k];
			*(Tuwpai+k) += Tuwpai2d[j][k];
			*(Tvwpai+k) += Tvwpai2d[j][k];
			
			*(Paiuu+k) += Paiuu2d[j][k];
			*(Paivv+k) += Paivv2d[j][k];
			*(Paiww+k) += Paiww2d[j][k];
			*(Paiuv+k) += Paiuv2d[j][k];
			*(Paiuw+k) += Paiuw2d[j][k];
			*(Paivw+k) += Paivw2d[j][k]; 
			
			*(Tiis+k) += Tiis2d[j][k];
			*(Tuus+k) += Tuus2d[j][k];
			*(Tvvs+k) += Tvvs2d[j][k];
			*(Twws+k) += Twws2d[j][k];
			*(Tuvs+k) += Tuvs2d[j][k];
			*(Tuws+k) += Tuws2d[j][k];
			*(Tvws+k) += Tvws2d[j][k];
			
			*(Tiir+k) += Tiir2d[j][k];
			*(Tuur+k) += Tuur2d[j][k];
			*(Tvvr+k) += Tvvr2d[j][k];
			*(Twwr+k) += Twwr2d[j][k];
			*(Tuvr+k) += Tuvr2d[j][k];
			*(Tuwr+k) += Tuwr2d[j][k];
			*(Tvwr+k) += Tvwr2d[j][k];
			
			*(Epii+k) += Epii2d[j][k];
			*(Epuu+k) += Epuu2d[j][k];
			*(Epvv+k) += Epvv2d[j][k];
			*(Epww+k) += Epww2d[j][k];
			*(Epuv+k) += Epuv2d[j][k];
			*(Epuw+k) += Epuw2d[j][k];
			*(Epvw+k) += Epvw2d[j][k];
			
			*(VdKdy+k) += VdKdy2d[j][k];
			*(WdKdz+k) += WdKdz2d[j][k];
			
			*(Vduudy+k) += Vduudy2d[j][k];
			*(Wduudz+k) += Wduudz2d[j][k];
			*(Vdvvdy+k) += Vdvvdy2d[j][k];
			*(Wdvvdz+k) += Wdvvdz2d[j][k];
			*(Vdwwdy+k) += Vdwwdy2d[j][k];
			*(Wdwwdz+k) += Wdwwdz2d[j][k];
			*(Vduvdy+k) += Vduvdy2d[j][k];
			*(Wduvdz+k) += Wduvdz2d[j][k];
			*(Vduwdy+k) += Vduwdy2d[j][k];
			*(Wduwdz+k) += Wduwdz2d[j][k];
			*(Vdvwdy+k) += Vdvwdy2d[j][k];
			*(Wdvwdz+k) += Wdvwdz2d[j][k];

			p1[k] += P1[j][k];
			p2[k] += P2[j][k];
			p3[k] += P3[j][k];
			p4[k] += P4[j][k];
			p5[k] += P5[j][k];
			p6[k] += P6[j][k];
			p7[k] += P7[j][k];
			p8[k] += P8[j][k];
			p9[k] += P9[j][k];
			p10[k] += P10[j][k];
			p11[k] += P11[j][k];
			p12[k] += P12[j][k];
			p13[k] += P13[j][k];
			p14[k] += P14[j][k];
			p15[k] += P15[j][k];
			p16[k] += P16[j][k];
			p17[k] += P17[j][k];
			p18[k] += P18[j][k];
		}
	}	
	
	for (k=0;k<zd;k++)
	{
		*(KE+k) /= ((double)exty);
	
		*(Pii+k) /= ((double)exty);
		*(Puu+k) /= ((double)exty);
		*(Pvv+k) /= ((double)exty);
		*(Pww+k) /= ((double)exty);
		*(Puv+k) /= ((double)exty);
		*(Puw+k) /= ((double)exty);
		*(Pvw+k) /= ((double)exty);
		
		*(Tiipai+k) /= ((double)exty);
		*(Tuupai+k) /= ((double)exty);
		*(Tvvpai+k) /= ((double)exty);
		*(Twwpai+k) /= ((double)exty);
		*(Tuvpai+k) /= ((double)exty);
		*(Tuwpai+k) /= ((double)exty);
		*(Tvwpai+k) /= ((double)exty);
		
		*(Paiuu+k) /= ((double)exty);
		*(Paivv+k) /= ((double)exty);
		*(Paiww+k) /= ((double)exty);
		*(Paiuv+k) /= ((double)exty);
		*(Paiuw+k) /= ((double)exty);
		*(Paivw+k) /= ((double)exty);		
		
		*(Tiis+k) /= ((double)exty);
		*(Tuus+k) /= ((double)exty);
		*(Tvvs+k) /= ((double)exty);
		*(Twws+k) /= ((double)exty);
		*(Tuvs+k) /= ((double)exty);
		*(Tuws+k) /= ((double)exty);
		*(Tvws+k) /= ((double)exty);
		
		*(Tiir+k) /= ((double)exty);
		*(Tuur+k) /= ((double)exty);
		*(Tvvr+k) /= ((double)exty);
		*(Twwr+k) /= ((double)exty);
		*(Tuvr+k) /= ((double)exty);
		*(Tuwr+k) /= ((double)exty);
		*(Tvwr+k) /= ((double)exty);
		
		*(Epii+k) /= ((double)exty);
		*(Epuu+k) /= ((double)exty);
		*(Epvv+k) /= ((double)exty);
		*(Epww+k) /= ((double)exty);
		*(Epuv+k) /= ((double)exty);
		*(Epuw+k) /= ((double)exty);
		*(Epvw+k) /= ((double)exty);
		
		*(VdKdy+k) /= ((double)exty);
		*(WdKdz+k) /= ((double)exty);
		
		*(Vduudy+k) /= ((double)exty);
		*(Wduudz+k) /= ((double)exty);
		*(Vdvvdy+k) /= ((double)exty);
		*(Wdvvdz+k) /= ((double)exty);
		*(Vdwwdy+k) /= ((double)exty);
		*(Wdwwdz+k) /= ((double)exty);
		*(Vduvdy+k) /= ((double)exty);
		*(Wduvdz+k) /= ((double)exty);
		*(Vduwdy+k) /= ((double)exty);
		*(Wduwdz+k) /= ((double)exty);
		*(Vdvwdy+k) /= ((double)exty);
		*(Wdvwdz+k) /= ((double)exty);

		p1[k] /= ((double)exty);
		p2[k] /= ((double)exty);
		p3[k] /= ((double)exty);
		p4[k] /= ((double)exty);
		p5[k] /= ((double)exty);
		p6[k] /= ((double)exty);
		p7[k] /= ((double)exty);
		p8[k] /= ((double)exty);
		p9[k] /= ((double)exty);
		p10[k] /= ((double)exty);
		p11[k] /= ((double)exty);
		p12[k] /= ((double)exty);
		p13[k] /= ((double)exty);
		p14[k] /= ((double)exty);
		p15[k] /= ((double)exty);
		p16[k] /= ((double)exty);
		p17[k] /= ((double)exty);
		p18[k] /= ((double)exty);
	}
	
	if (pd.myrank==0)
	{
		sprintf(fn,"turb-fld-2d.%.3d.%.4d",pd.myrank,(int)ts);
		sv=fopen(fn,"wb");
		fwrite(UAVG,sizeof(DP),zd,sv);
		fwrite(VAVG,sizeof(DP),zd,sv);
		fwrite(WAVG,sizeof(DP),zd,sv);
		fwrite(UTIN,sizeof(DP),zd,sv);
		fwrite(VTIN,sizeof(DP),zd,sv);
		fwrite(WTIN,sizeof(DP),zd,sv);
		fwrite(uv,sizeof(DP),zd,sv);
		fwrite(uw,sizeof(DP),zd,sv);
		fwrite(vw,sizeof(DP),zd,sv);
		fwrite(DAVG,sizeof(DP),zd,sv);
		fwrite(DTIN,sizeof(DP),zd,sv);
		fwrite(SKEWu,sizeof(DP),zd,sv);
		fwrite(SKEWv,sizeof(DP),zd,sv);
		fwrite(SKEWw,sizeof(DP),zd,sv);
		fwrite(SKEWuv,sizeof(DP),zd,sv);
		fwrite(SKEWuw,sizeof(DP),zd,sv);
		fwrite(SKEWvw,sizeof(DP),zd,sv);
		fwrite(SKEWp,sizeof(DP),zd,sv);
		fwrite(KURTu,sizeof(DP),zd,sv);
		fwrite(KURTv,sizeof(DP),zd,sv);
		fwrite(KURTw,sizeof(DP),zd,sv);
		fwrite(KURTuv,sizeof(DP),zd,sv);
		fwrite(KURTuw,sizeof(DP),zd,sv);
		fwrite(KURTvw,sizeof(DP),zd,sv);
		fwrite(KURTp,sizeof(DP),zd,sv);
		fwrite(UW,sizeof(DP),zd,sv);
		fclose(sv);

		sprintf(fn,"vel-field-2d.%.3d.%.4d",pd.myrank,(int)ts);
		sv=fopen(fn,"wb");
		fwrite(uu2d,sizeof(DP),(yd*zd),sv);
		fwrite(vv2d,sizeof(DP),(yd*zd),sv);
		fwrite(ww2d,sizeof(DP),(yd*zd),sv);
		fwrite(uv2d,sizeof(DP),(yd*zd),sv);
		fwrite(uw2d,sizeof(DP),(yd*zd),sv);
		fwrite(vw2d,sizeof(DP),(yd*zd),sv);
		fwrite(UW2d,sizeof(DP),(yd*zd),sv);
		fwrite(dd2d,sizeof(DP),(yd*zd),sv);
		fwrite(dAVG,sizeof(DP),(yd*zd),sv);
		fclose(sv);
		
		sprintf(fn,"vort-fluc-2d.%.3d.%.4d",pd.myrank,(int)ts);	
		sv=fopen(fn,"wb");
		fwrite(wx,sizeof(DP),zd,sv);
		fwrite(wy,sizeof(DP),zd,sv);
		fwrite(wz,sizeof(DP),zd,sv);
		fclose(sv);
		
		sprintf(fn,"vort-fluc-2d-field.%.3d.%.4d",pd.myrank,(int)ts);	
		sv=fopen(fn,"wb");
		fwrite(&omegax2d[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&omegay2d[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&omegaz2d[0][0],sizeof(DP),(zd*yd),sv);
		fclose(sv);
		
		sprintf(fn,"Mean-X-Vorticity-2d.%.3d.%.4d",pd.myrank,(int)ts);	
		sv=fopen(fn,"wb");
		fwrite(&OmegaX[0][0],sizeof(DP),yd*zd,sv);
		fclose(sv);
		
		sprintf(fn,"uAVG-2d.%.3d.%.4d",pd.myrank,(int)ts);
		sv=fopen(fn,"wb");
		fwrite(&uAVG[0][0],sizeof(DP),yd*zd,sv);
		fclose(sv);

		sprintf(fn,"vAVG-2d.%.3d.%.4d",pd.myrank,(int)ts);
        sv=fopen(fn,"wb");
        fwrite(&vAVG[0][0],sizeof(DP),yd*zd,sv);
        fclose(sv);

		sprintf(fn,"wAVG-2d.%.3d.%.4d",pd.myrank,(int)ts);
        sv=fopen(fn,"wb");
        fwrite(&wAVG[0][0],sizeof(DP),yd*zd,sv);
        fclose(sv);
                
        sprintf(fn,"KE-Budget-2d.%.3d.%.4d",pd.myrank,(int)ts);
        sv=fopen(fn,"wb");
		fwrite(Pii,sizeof(DP),zd,sv);
		fwrite(Puu,sizeof(DP),zd,sv);
		fwrite(Pvv,sizeof(DP),zd,sv);
		fwrite(Pww,sizeof(DP),zd,sv);
		fwrite(Puv,sizeof(DP),zd,sv);
		fwrite(Puw,sizeof(DP),zd,sv);
		fwrite(Pvw,sizeof(DP),zd,sv);
		fwrite(Tiir,sizeof(DP),zd,sv);
		fwrite(Tuur,sizeof(DP),zd,sv);
		fwrite(Tvvr,sizeof(DP),zd,sv);
		fwrite(Twwr,sizeof(DP),zd,sv);
		fwrite(Tuvr,sizeof(DP),zd,sv);
		fwrite(Tuwr,sizeof(DP),zd,sv);
		fwrite(Tvwr,sizeof(DP),zd,sv);
		fwrite(Tiipai,sizeof(DP),zd,sv);
		fwrite(Tuupai,sizeof(DP),zd,sv);
		fwrite(Tvvpai,sizeof(DP),zd,sv);
		fwrite(Twwpai,sizeof(DP),zd,sv);
		fwrite(Tuvpai,sizeof(DP),zd,sv);
		fwrite(Tuwpai,sizeof(DP),zd,sv);
		fwrite(Tvwpai,sizeof(DP),zd,sv);
		fwrite(Paiuu,sizeof(DP),zd,sv);
		fwrite(Paivv,sizeof(DP),zd,sv);
		fwrite(Paiww,sizeof(DP),zd,sv);
		fwrite(Paiuv,sizeof(DP),zd,sv);
		fwrite(Paiuw,sizeof(DP),zd,sv);
		fwrite(Paivw,sizeof(DP),zd,sv);
		fwrite(Tiis,sizeof(DP),zd,sv);
		fwrite(Tuus,sizeof(DP),zd,sv);
		fwrite(Tvvs,sizeof(DP),zd,sv);
		fwrite(Twws,sizeof(DP),zd,sv);
		fwrite(Tuvs,sizeof(DP),zd,sv);
		fwrite(Tuws,sizeof(DP),zd,sv);
		fwrite(Tvws,sizeof(DP),zd,sv);
		fwrite(Epii,sizeof(DP),zd,sv);
		fwrite(Epuu,sizeof(DP),zd,sv);
		fwrite(Epvv,sizeof(DP),zd,sv);
		fwrite(Epww,sizeof(DP),zd,sv);
		fwrite(Epuv,sizeof(DP),zd,sv);
		fwrite(Epuw,sizeof(DP),zd,sv);
		fwrite(Epvw,sizeof(DP),zd,sv);
		fwrite(KE,sizeof(DP),zd,sv);
		fwrite(VdKdy,sizeof(DP),zd,sv);
		fwrite(WdKdz,sizeof(DP),zd,sv);
		fwrite(Vduudy,sizeof(DP),zd,sv);
		fwrite(Wduudz,sizeof(DP),zd,sv);
		fwrite(Vdvvdy,sizeof(DP),zd,sv);
		fwrite(Wdvvdz,sizeof(DP),zd,sv);
		fwrite(Vdwwdy,sizeof(DP),zd,sv);
		fwrite(Wdwwdz,sizeof(DP),zd,sv);
		fwrite(Vduvdy,sizeof(DP),zd,sv);
		fwrite(Wduvdz,sizeof(DP),zd,sv);
		fwrite(Vduwdy,sizeof(DP),zd,sv);
		fwrite(Wduwdz,sizeof(DP),zd,sv);
		fwrite(Vdvwdy,sizeof(DP),zd,sv);
		fwrite(Wdvwdz,sizeof(DP),zd,sv);
		fclose(sv);

		sprintf(fn,"TP-Breakdown.%.3d.%.4d",pd.myrank,(int)ts);
		sv=fopen(fn,"wb");
		fwrite(&p1[0],sizeof(DP),zd,sv);
		fwrite(&p2[0],sizeof(DP),zd,sv);
		fwrite(&p3[0],sizeof(DP),zd,sv);
		fwrite(&p4[0],sizeof(DP),zd,sv);
		fwrite(&p5[0],sizeof(DP),zd,sv);
		fwrite(&p6[0],sizeof(DP),zd,sv);
		fwrite(&p7[0],sizeof(DP),zd,sv);
		fwrite(&p8[0],sizeof(DP),zd,sv);
		fwrite(&p9[0],sizeof(DP),zd,sv);
		fwrite(&p10[0],sizeof(DP),zd,sv);
		fwrite(&p11[0],sizeof(DP),zd,sv);
		fwrite(&p12[0],sizeof(DP),zd,sv);
		fwrite(&p13[0],sizeof(DP),zd,sv);
		fwrite(&p14[0],sizeof(DP),zd,sv);
		fwrite(&p15[0],sizeof(DP),zd,sv);
		fwrite(&p16[0],sizeof(DP),zd,sv);
		fwrite(&p17[0],sizeof(DP),zd,sv);
		fwrite(&p18[0],sizeof(DP),zd,sv);
		fclose(sv);

		sprintf(fn,"TP-Breakdown-2d-field.%.3d.%.4d",pd.myrank,(int)ts);
		sv=fopen(fn,"wb");
		fwrite(&P1[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P2[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P3[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P4[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P5[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P6[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P7[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P8[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P9[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P10[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P11[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P12[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P13[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P14[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P15[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P16[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P17[0][0],sizeof(DP),(zd*yd),sv);
		fwrite(&P18[0][0],sizeof(DP),(zd*yd),sv);
		fclose(sv);


		sprintf(fn,"Budget-2d-field.%.3d.%.4d",pd.myrank,(int)ts);
		sv=fopen(fn,"wb");
		fwrite(Pii2d,sizeof(DP),(zd*yd),sv);
		fwrite(Puu2d,sizeof(DP),(zd*yd),sv);
		fwrite(Pvv2d,sizeof(DP),(zd*yd),sv);
		fwrite(Pww2d,sizeof(DP),(zd*yd),sv);
		fwrite(Puv2d,sizeof(DP),(zd*yd),sv);
		fwrite(Puw2d,sizeof(DP),(zd*yd),sv);
		fwrite(Pvw2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tiir2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tuur2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tvvr2d,sizeof(DP),(zd*yd),sv);
		fwrite(Twwr2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tuvr2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tuwr2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tvwr2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tiipai2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tuupai2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tvvpai2d,sizeof(DP),(zd*yd),sv);
		fwrite(Twwpai2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tuvpai2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tuwpai2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tvwpai2d,sizeof(DP),(zd*yd),sv);
		fwrite(Paiuu2d,sizeof(DP),(zd*yd),sv);
		fwrite(Paivv2d,sizeof(DP),(zd*yd),sv);
		fwrite(Paiww2d,sizeof(DP),(zd*yd),sv);
		fwrite(Paiuv2d,sizeof(DP),(zd*yd),sv);
		fwrite(Paiuw2d,sizeof(DP),(zd*yd),sv);
		fwrite(Paivw2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tiis2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tuus2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tvvs2d,sizeof(DP),(zd*yd),sv);
		fwrite(Twws2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tuvs2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tuws2d,sizeof(DP),(zd*yd),sv);
		fwrite(Tvws2d,sizeof(DP),(zd*yd),sv);
		fwrite(Epii2d,sizeof(DP),(zd*yd),sv);
		fwrite(Epuu2d,sizeof(DP),(zd*yd),sv);
		fwrite(Epvv2d,sizeof(DP),(zd*yd),sv);
		fwrite(Epww2d,sizeof(DP),(zd*yd),sv);
		fwrite(Epuv2d,sizeof(DP),(zd*yd),sv);
		fwrite(Epuw2d,sizeof(DP),(zd*yd),sv);
		fwrite(Epvw2d,sizeof(DP),(zd*yd),sv);
		fwrite(KE2d,sizeof(DP),(zd*yd),sv);
		fwrite(VdKdy2d,sizeof(DP),(zd*yd),sv);
		fwrite(WdKdz2d,sizeof(DP),(zd*yd),sv);
		fwrite(Vduudy2d,sizeof(DP),(zd*yd),sv);
		fwrite(Wduudz2d,sizeof(DP),(zd*yd),sv);
		fwrite(Vdvvdy2d,sizeof(DP),(zd*yd),sv);
		fwrite(Wdvvdz2d,sizeof(DP),(zd*yd),sv);
		fwrite(Vdwwdy2d,sizeof(DP),(zd*yd),sv);
		fwrite(Wdwwdz2d,sizeof(DP),(zd*yd),sv);
		fwrite(Vduvdy2d,sizeof(DP),(zd*yd),sv);
		fwrite(Wduvdz2d,sizeof(DP),(zd*yd),sv);
		fwrite(Vduwdy2d,sizeof(DP),(zd*yd),sv);
		fwrite(Wduwdz2d,sizeof(DP),(zd*yd),sv);
		fwrite(Vdvwdy2d,sizeof(DP),(zd*yd),sv);
		fwrite(Wdvwdz2d,sizeof(DP),(zd*yd),sv);
		fclose(sv);
		
		sprintf(fn,"dudzusuw-2d.%.3d.%.4d",pd.myrank,(int)ts);
                sv=fopen(fn,"wb");
                fwrite(&dudzusuw[0][0],sizeof(DP),yd*zd,sv);
                fclose(sv);
		
		sprintf(fn,"dquadrant.%.3d.%.4d",pd.myrank,(int)ts);
		sv=fopen(fn,"wb");
		fwrite(&quadrant[0][0],sizeof(DP),4*zd,sv);
		fwrite(uw,sizeof(DP),zd,sv);
		fclose(sv);
		
		sprintf(fn,"dquadrant-2d.%.3d.%.4d",pd.myrank,(int)ts);
		sv=fopen(fn,"wb");
		fwrite(&quadrant2d[0][0][0],sizeof(DP),4*yd*zd,sv);
		fwrite(uw,sizeof(DP),zd,sv);
		fwrite(&uw2d[0][0],sizeof(DP),yd*zd,sv);
		fclose(sv);

	}
	
	
	free(UAVG);
	free(VAVG);
	free(WAVG);
	free(DAVG);
	free(DUMMY);
	free(UW);
	
	free(UTIN);
	free(VTIN);
	free(WTIN);
	free(DTIN);
	
	free(DUMMYYZ);
	
	free(uu);
	free(vv);
	free(ww);
	free(uw);
	free(vw);
	free(uv);
	free(dd);

	free(Pii);
	free(Puu);
	free(Pvv);
	free(Pww);
	free(Puw);
	free(Puv);
	free(Pvw);
	
	free(Tiir);
	free(Tiipai);
	free(Tiis);
	free(Epii);
	
 	free(KE);

	free(Paiuu);
	free(Paivv);
	free(Paiww);
	free(Paiuw);
	free(Paiuv);
	free(Paivw);

	free(Tiirt);
	free(Tuur);
	free(Tvvr);
	free(Twwr);
	free(Tuwr);
	free(Tuvr);
	free(Tvwr);

	free(Tuus);
	free(Tvvs);
	free(Twws);
	free(Tuws);
	free(Tuvs);
	free(Tvws);

	free(Tuupai);
	free(Tvvpai);
	free(Twwpai);
	free(Tuwpai);
	free(Tuvpai);
	free(Tvwpai);

	free(Epuu);
	free(Epvv);
	free(Epww);
	free(Epuw);
	free(Epuv);
	free(Epvw);
	
	free(wx);
	free(wy);
	free(wz);
	
	free(SKEWu);
	free(SKEWv);
	free(SKEWw);
	free(SKEWp);
	free(SKEWuw);
	free(SKEWuv);
	free(SKEWvw);
	
	free(KURTu);
	free(KURTv);
	free(KURTw);
	free(KURTp);
	free(KURTuw);
	free(KURTuv);
	free(KURTvw);
	
	free(VdKdy);
	free(WdKdz);
	
	free(Vduudy);
	free(Wduudz);
	free(Vdvvdy);
	free(Wdvvdz);
	free(Vdwwdy);
	free(Wdwwdz);
	free(Vduvdy);
	free(Wduvdz);
	free(Vduwdy);
	free(Wduwdz);
	free(Vdvwdy);
	free(Wdvwdz);
	
	free(iDUMMYq);
	free(DUMMYq2d);
	free(DUMMYq);
	
return V;
}