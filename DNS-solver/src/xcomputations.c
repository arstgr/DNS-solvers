#include "definitions.h"

/*Performing Collision and streaming Parallel cut planes normal to x direction */

POINTER xcomputations(int xd, int yd, int zd, POINTER V,PDATA pd,DP Gx, DP Gy, DP Gz, DP rhozero,DP tau, DP ubulk)
{
	int i,j,k,a,t,x[4]={0,1,(xd-2),(xd-1)},y[4]={0,1,(yd-2),(yd-1)};
	int ito,jto,kto;
	int zt[4]={0,1,(zd-2),(zd-1)},num=3,ct;
	const int L[5]={2,8,9,13,14};
	const int R[5]={1,7,10,11,12};
	const int U[5]={3,7,8,15,16};
	const int D[5]={4,9,10,17,18};
	const int UR[9]={1,3,5,6,7,11,12,15,16};
	const int DL[9]={2,4,5,6,9,13,14,17,18};
	const int UL[9]={2,3,5,6,8,13,14,15,16};
	const int DR[9]={1,4,5,6,10,11,12,17,18};
	int tag,flag;
	MPI_Status status[12]={0},tstatus;
	DP v;
	DP feq,u,z,dens,*ptr,*adr;
	int solid,sld1,sld2,sld;
	DP ux,uy,uz,ueq[3];
	int ii,jj,kk,il,jl,kl,ic,jc,kc,xc,yc,zc;
	DP Fi,*taddrs;
	DP rl,rr,dq;
	
	DP tpsend[zd][9],tprecv[zd][9],tprecv2[zd][9];

	Gx=V.force;
	for (t=0;t<4;t++)
	{
		i=x[t];
		for (j=0;j<yd;j++)
		{
			for (k=0;k<zd;k++)
			{
				solid = Fs(V.s,i,j,k,xd,yd,zd);
				if (solid>-1)
				{	
					v = 0.;
					ux=0.;uy=0.;uz=0.;dens=0.;
					for (a=0;a<19;a++)
					{
						ux += E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						uy += E1[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						uz += E2[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						dens += Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					}
					ueq[0]=ux/dens;
					ueq[1]=uy/dens;
					ueq[2]=uz/dens;
					v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
					
					for (a=0;a<19;a++)
					{
						u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
						feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
						//Fi=w[a]*(1-0.5*tau)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						//Fi=BI[a]*dens*(E0[a]*Gx+E1[a]*Gy+E2[a]*Gz);
						Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						//Fi=E0[a]*Gx/12.;
						//Fi=BI[a]*E0[a]*Gx;
						//Fi=0.;
						Fb(V.ftemp,i,j,k,a,xd,yd,zd,19) += (Fi-tau*(Fb(V.ftemp,i,j,k,a,xd,yd,zd,19)-feq));
					}
				}
			}
		}
	}
	for (i=2;i<(xd-2);i++)
	{
		for (t=0;t<4;t++)
		{
			j=y[t];
			for (k=0;k<zd;k++)
			{
				solid = Fs(V.s,i,j,k,xd,yd,zd);
				if (solid>-1)
				{	
					v = 0.;
					ux=0.;uy=0.;uz=0.;dens=0.;
					for (a=0;a<19;a++)
					{
						ux += E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						uy += E1[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						uz += E2[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						dens += Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					}
					ueq[0]=ux/dens;
					ueq[1]=uy/dens;
					ueq[2]=uz/dens;
					v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
					
					for (a=0;a<19;a++)
					{
						u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
						feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
						//Fi=w[a]*(1-0.5*tau)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						//Fi=BI[a]*dens*(E0[a]*Gx+E1[a]*Gy+E2[a]*Gz);
						Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
						//Fi=E0[a]*Gx/12.;
						//Fi=BI[a]*E0[a]*Gx;
						//Fi=0.;
						Fb(V.ftemp,i,j,k,a,xd,yd,zd,19) += (Fi-tau*(Fb(V.ftemp,i,j,k,a,xd,yd,zd,19)-feq));
					}
				}
			}
		}
	}
	for (t=0;t<4;t++)
	{
		i=x[t];
		for (j=0;j<yd;j++)
		{
			for (k=0;k<zd;k++)
			{
				solid = Fs(V.s,i,j,k,xd,yd,zd);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						sld=Fs(V.s,ito,jto,kto,xd,yd,zd);
						sld1=Fs(V.s,i,j,kto,xd,yd,zd);
						if (sld>-1)
							Fb(V.f,ito,jto,kto,a,xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						else if ((sld==-1)&&(sld1==-1))
							Fb(V.f,i,j,k,opposite[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						else if ((sld==-1)&&(sld1==-2))
							Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						else if (sld==-2)
							Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((ito>-1)&&(ito<xd)&&(jto>-1)&&(jto<yd)&&(kto>-1)&&(kto<zd))
						{
							sld=Fs(V.s,ito,jto,kto,xd,yd,zd);
							sld1=Fs(V.s,i,j,kto,xd,yd,zd);
							if (sld>-1)
								Fb(V.f,ito,jto,kto,a,xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
							else if ((sld==-1)&&(sld1==-1))
								Fb(V.f,i,j,k,opposite[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
							else if ((sld==-1)&&(sld1==-2))
								Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
							else if (sld==-2)
								Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						}
					}
				}
			}
		}
	}
	for (i=2;i<(xd-2);i++)
	{
		for (t=0;t<4;t++)
		{
			j=y[t];
			for (k=0;k<zd;k++)
			{
				solid = Fs(V.s,i,j,k,xd,yd,zd);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						sld=Fs(V.s,ito,jto,kto,xd,yd,zd);
						sld1=Fs(V.s,i,j,kto,xd,yd,zd);
						if (sld>-1)
							Fb(V.f,ito,jto,kto,a,xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						else if ((sld==-1)&&(sld1==-1))
							Fb(V.f,i,j,k,opposite[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						else if ((sld==-1)&&(sld1==-2))
							Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						else if (sld==-2)
							Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((ito>-1)&&(ito<xd)&&(jto>-1)&&(jto<yd)&&(kto>-1)&&(kto<zd))
						{
							sld=Fs(V.s,ito,jto,kto,xd,yd,zd);
							sld1=Fs(V.s,i,j,kto,xd,yd,zd);
							if (sld>-1)
								Fb(V.f,ito,jto,kto,a,xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
							else if ((sld==-1)&&(sld1==-1))
								Fb(V.f,i,j,k,opposite[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
							else if ((sld==-1)&&(sld1==-2))
								Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
							else if (sld==-2)
								Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						}
					}
				}
			}
		}
	}

	i=0;
	ptr=V.leftbufs;
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			for (a=0;a<5;a++)
			{
				*(ptr++)=Fb(V.f,i,j,k,L[a],xd,yd,zd,19);
			}
		}
	}

	i=xd-1;
	adr=V.rightbufs;
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			for (a=0;a<5;a++)
			{
				(*(adr++))=Fb(V.f,i,j,k,R[a],xd,yd,zd,19);
			}
		}
	}
	
	j=0;
	ptr=V.dwbufs;
	for (i=0;i<xd;i++)
	{
		for (k=0;k<zd;k++)
		{
			for (a=0;a<5;a++)
			{
				*(ptr++)=Fb(V.f,i,j,k,D[a],xd,yd,zd,19);
			}
		}
	}

	j=yd-1;
	adr=V.upbufs;
	for (i=0;i<xd;i++)
	{
		for (k=0;k<zd;k++)
		{
			for (a=0;a<5;a++)
			{
				(*(adr++))=Fb(V.f,i,j,k,U[a],xd,yd,zd,19);
			}
		}
	}
	
	i=0;j=0;
	adr=V.dlbufs;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			(*(adr++))=Fb(V.f,i,j,k,DL[a],xd,yd,zd,19);
	}
	i=xd-1;j=yd-1;
	adr=V.urbufs;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			(*(adr++))=Fb(V.f,i,j,k,UR[a],xd,yd,zd,19);
	}
	
	MPI_Startall(12,&V.req[0]);
	
	il=(xd-4)/cs;
	jl=(yd-4)/cs;
	kl=zd/cs;

	xc=(xd-4)%cs;
	yc=(yd-4)%cs;
	zc=zd%cs;
	
	for (ii=0;ii<=il;ii++)
	{
		for (jj=0;jj<=jl;jj++)
		{
			for (kk=0;kk<=kl;kk++)
			{
				for (ic=0;(ic<cs)&&((i=ic+cs*ii+2)<(xd-2));ic++)
				{
					for (jc=0;(jc<cs)&&((j=jc+cs*jj+2)<(yd-2));jc++)
					{
						for (kc=0;(kc<cs)&&((k=kc+cs*kk)<zd);kc++)
						{
							solid = Fs(V.s,i,j,k,xd,yd,zd);
							if (solid>-1)
							{	
								v = 0.;
								ux=0.;uy=0.;uz=0.;dens=0.;
								for (a=0;a<19;a++)
								{
									ux += E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
									uy += E1[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
									uz += E2[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
									dens += Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
								}
								ueq[0]=ux/dens;
								ueq[1]=uy/dens;
								ueq[2]=uz/dens;
								v = ueq[0]*ueq[0]+ueq[1]*ueq[1]+ueq[2]*ueq[2];
								
								for (a=0;a<19;a++)
								{
									u=E0[a]*ueq[0]+E1[a]*ueq[1]+E2[a]*ueq[2];
									feq=w[a]*dens*(1.+3.*u+4.5*u*u-1.5*v);
									//Fi=w[a]*(1-0.5*tau)*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
									//Fi=BI[a]*dens*(E0[a]*Gx+E1[a]*Gy+E2[a]*Gz);
									//Fi=3.*dens*w[a]*E0[a]*Gx;
									//Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
									Fi=w[a]*(3.*(E0[a]-ueq[0])+9.*u*E0[a])*Gx;
									//Fi=3.*dens*(1-0.5*tau)*w[a]*E0[a]*Gx;
									//Fi=E0[a]*Gx/12.;
									//Fi=BI[a]*E0[a]*Gx;
									//Fi=0.;
									Fb(V.ftemp,i,j,k,a,xd,yd,zd,19) += (Fi-tau*(Fb(V.ftemp,i,j,k,a,xd,yd,zd,19)-feq));
								}
							}
						}
					}
				}
			}
		}
	}
	MPI_Testall(12,&V.req[0],&flag,&status[0]);
	for (i=2;i<(xd-2);i++)
	{
		for (j=2;j<(yd-2);j++)
		{
			for (t=0;t<4;t++)
			{
				k=zt[t];
				solid = Fs(V.s,i,j,k,xd,yd,zd);
				if (solid==0)
				{
					for (a=0;a<19;a++)
					{
						ito=i+e0[a];
						jto=j+e1[a];
						kto=k+e2[a];
						sld=Fs(V.s,ito,jto,kto,xd,yd,zd);
						sld1=Fs(V.s,i,j,kto,xd,yd,zd);
						if (sld>-1)
							Fb(V.f,ito,jto,kto,a,xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						else if ((sld==-1)&&(sld1==-1))
							Fb(V.f,i,j,k,opposite[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						else if ((sld==-1)&&(sld1==-2))
							Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						else if (sld==-2)
							Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
					}
				}
				else if (solid==1)
				{
					for (a=0;a<19;a++)
					{
						jto=j+e1[a];
						
						kto=k+e2[a];
						
						ito=i+e0[a];
						
						if ((kto>-1)&&(kto<zd))
						{
							sld=Fs(V.s,ito,jto,kto,xd,yd,zd);
							sld1=Fs(V.s,i,j,kto,xd,yd,zd);
							if (sld>-1)
								Fb(V.f,ito,jto,kto,a,xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
							else if ((sld==-1)&&(sld1==-1))
								Fb(V.f,i,j,k,opposite[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
							else if ((sld==-1)&&(sld1==-2))
								Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
							else if (sld==-2)
								Fb(V.f,ito,jto,k,mirror[a],xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
						}
					}
				}
			}
		}
	}
	for (ii=0;ii<=il;ii++)
	{
		for (jj=0;jj<=jl;jj++)
		{
			for (kk=0;kk<=kl;kk++)
			{
				for (ic=0;(ic<cs)&&((i=ic+cs*ii+2)<(xd-2));ic++)
				{
					for (jc=0;(jc<cs)&&((j=jc+cs*jj+2)<(yd-2));jc++)
					{
						for (kc=0;(kc<cs)&&((k=kc+cs*kk+2)<(zd-2));kc++)
						{
							for (a=0;a<19;a++)
							{
								ito=i+e0[a];
								jto=j+e1[a];

								kto=k+e2[a];
								
								Fb(V.f,ito,jto,kto,a,xd,yd,zd,19)=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
							}
						}
					}
				}
			}
			MPI_Testall(12,&V.req[0],&flag,&status[0]);
		}
	}

	MPI_Waitall(12,&V.req[0],&status[0]);

	i=0;
	adr=V.leftbufr;
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.f,i,j,k,R[a],xd,yd,zd,19)=(*(adr++));
			}
		}
	}

	ptr=V.rightbufr;
	i=xd-1;
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.f,i,j,k,L[a],xd,yd,zd,19)=(*(ptr++));
			}
		}
	}
	
	j=0;
	adr=V.dwbufr;
	for (i=0;i<xd;i++)
	{
		for (k=0;k<zd;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.f,i,j,k,U[a],xd,yd,zd,19)=(*(adr++));
			}
		}
	}

	ptr=V.upbufr;
	j=yd-1;
	for (i=0;i<xd;i++)
	{
		for (k=0;k<zd;k++)
		{
			for (a=0;a<5;a++)
			{
				Fb(V.f,i,j,k,D[a],xd,yd,zd,19)=(*(ptr++));
			}
		}
	}
	
	i=0;j=0;
	ptr=V.dlbufr;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,UR[a],xd,yd,zd,19)=(*(ptr++));
	}

	i=xd-1;j=yd-1;
	ptr=V.urbufr;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,DL[a],xd,yd,zd,19)=(*(ptr++));
	}
	
	i=0;j=0;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.f,i,j,k,UR[a],xd,yd,zd,19);
	}
	

	MPI_Sendrecv(tpsend,9*zd,MPI_DOUBLE,pd.left,221,tprecv,9*zd,MPI_DOUBLE,pd.right,221,MPI_COMM_WORLD,&tstatus);
	MPI_Sendrecv(tpsend,9*zd,MPI_DOUBLE,pd.dw,222,tprecv2,9*zd,MPI_DOUBLE,pd.up,222,MPI_COMM_WORLD,&tstatus);
	
	i=xd-1;j=0;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,UR[a],xd,yd,zd,19) = tprecv[k][a];
	}
	i=0;j=yd-1;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,UR[a],xd,yd,zd,19) = tprecv2[k][a];
	}
	
	i=0;j=yd-1;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.f,i,j,k,UL[a],xd,yd,zd,19);
	}
	MPI_Sendrecv(tpsend,9*zd,MPI_DOUBLE,pd.ul,223,tprecv,9*zd,MPI_DOUBLE,pd.dr,223,MPI_COMM_WORLD,&tstatus);
	i=xd-1;j=0;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,UL[a],xd,yd,zd,19) = tprecv[k][a];
	}
	
	i=xd-1;j=0;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			tpsend[k][a] = Fb(V.f,i,j,k,DR[a],xd,yd,zd,19);
	}
	MPI_Sendrecv(tpsend,9*zd,MPI_DOUBLE,pd.dr,224,tprecv,9*zd,MPI_DOUBLE,pd.ul,224,MPI_COMM_WORLD,&tstatus);
	i=0;j=yd-1;
	for (k=0;k<zd;k++)
	{
		for (a=0;a<9;a++)
			Fb(V.f,i,j,k,DR[a],xd,yd,zd,19) = tprecv[k][a];
	}
	
	taddrs=V.f;
	V.f=V.ftemp;
	V.ftemp=taddrs;

return V;
}