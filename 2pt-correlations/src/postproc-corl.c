/********************************************************************************
 * Amirreza rastegari                                                           *
 * arstgri@gmail.com                                                            *
 * Creates the required files to plot 2 point velocity correlations             *
 * the output files, *.dat, can be read by general purpose plotting software    *
 * such as tecplot(R)                                                           *
 * It should be used with correlations.c                                        *
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define Reb 7200.0000
#define ubulk 0.0866034
#define TSTR 500
#define TEND 900
#define ZPTGT 170    /* Z+ for the requested plane */
//if defined has priority
//#define ZOH 0.8   /* Z/H for the requested plane */
#define DT 1
#define TEXT "1DAVG"
#define SLPTXT "1DAVG"

// width of the shear free area
//#define shfw 64

//width of slip
#define gplus 0
//width of no-slip
#define wplus 0
//half width of no-slip
#define gst 0

typedef double DP;

#define xdv 512
#define ydv 256
#define zdv 221

#define Fy(a,j,k,b,yd,zd,d) (*(a+(k*yd*d)+(j*d)+b))
#define Fx(a,i,k,b,xd,zd,d) (*(a+(k*xd*d)+(i*d)+b))

void correl_res(int xd, int yd, int zd, int st, int en)
{
	int i,j,k,a;
	FILE *sv;
	FILE *tcp,*wstr;
	char fn[100];
	int num=(en-st)+1;
	DP *z,dplus,Ret;
	DP utauavg=0.,utvt,time,nu;
	int cnt=0,ct=0,rd;
	int up;
	DP ubavg=0.,ubt=0.,twall,val,tau,ub=0.;
	char stri[60];
	DP uttemp,uttemp1,ttemp;
	double zn,zp;
	int nz;
	double *ycorl,*xcorl,*xcorlsl,*xcorlnsl,*ycorlt,*xcorlt,*xcorlslt,*xcorlnslt;
	double uavg[zd],uavgt[zd];
	int extx,exty;
	
	xcorl = (DP *)calloc((xdv*zd*3),sizeof(DP));
	xcorlsl = (DP *)calloc((xdv*zd*3),sizeof(DP));
	xcorlnsl = (DP *)calloc((xdv*zd*3),sizeof(DP));
	ycorl = (DP *)calloc((yd*zd*3),sizeof(DP));
	xcorlt = (DP *)calloc((xdv*zd*3),sizeof(DP));
	xcorlslt = (DP *)calloc((xdv*zd*3),sizeof(DP));
	xcorlnslt = (DP *)calloc((xdv*zd*3),sizeof(DP));
	ycorlt = (DP *)calloc((yd*zd*3),sizeof(DP));
	
	up=1+(zd-1)/2;
	extx = xdv/2;
	exty = yd/2;
	num/=DT;
	z=(DP *)calloc(zd,sizeof(DP));

	wstr=fopen("p-grad.txt","r");
	ttemp=TSTR;
	uttemp=0.;
	while (fscanf(wstr,"%lf %lf\n",&ttemp,&uttemp1)!=EOF)
	{
		if ((ttemp>=TSTR)&&(ttemp<=TEND))
		{
			uttemp+=uttemp1;
			ct++;
		}
	}
	fclose(wstr);

	utauavg=sqrt(uttemp/ct);
	tau=(3./Reb)*ubulk*(zd-2.)+0.5;
	nu=(1./6.)*(2.*tau-1.);
	Ret=utauavg*(zd-2.)*0.5/nu;
	dplus=Ret/(0.5*(zd-2));
	printf("ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);

	ct=0;ub=0.;
	sv=fopen("flowrate.txt","r");
	while (fscanf(sv,"%lf %lf\n",&ttemp,&ubt)!=EOF)
        {
                if ((ttemp>=TSTR)&&(ttemp<=TEND))
                {
                        ub+=ubt;
                        ct++;
                }
        }
	fclose(sv);
	ub=ub/((double)(ct*yd*(zd-2.)));
	printf("ub=%lf\n",ub);

	for (k=1;k<(zd-1);k++)
	{
		*(z+k)=(k-0.5*(zd-1))*dplus;
	}
	*z=(0.5-0.5*(zd-1))*dplus;
	*(z+zd-1)=(0.5*(zd-1)-0.5)*dplus;

	for (k=0;k<zd;k++)
	{
		*(z+k) /= Ret;
	}

	for (k=0;k<zd;k++)
		uavg[k]=0.;
	
	zn = -1;
	k=0;
	while (zn<ZPTGT)
	{
		zp = zn;
		zn = ((*(z+k))+1)*Ret;
		k++;
	}
	nz = k-2;
	if (fabs(zn-ZPTGT)<fabs(zp-ZPTGT))
	{
		zp = zn;
		nz++;
	}
	printf("nz=%d z+=%f\n",nz,zp); 

#ifdef ZOH
	zn=-1;
	k=0;
	while (zn<ZOH)
	{
		zp =zn;
		zn = ((*(z+k))+1);
		k++;
	}
	nz = k-2;
	if (fabs(zn-ZOH)<fabs(zp-ZOH))
	{
		zp = zn;
		nz++;
	}
	printf("forget about the previous one, nz=%d z/h=%f\n",nz,zp); 
#endif

	for (k=st;k<(en+1);k+=DT)
	{
		sprintf(fn,"xcorl.%.4d",k);
		sv=fopen(fn,"rb");
		fread(xcorlt,sizeof(DP),(zd*xdv*3),sv);
		fclose(sv);
/*
		sprintf(fn,"xcorlsl.%.4d",k);
		sv=fopen(fn,"rb");
		fread(xcorlslt,sizeof(DP),(zd*xdv*3),sv);
		fclose(sv);
		
		sprintf(fn,"xcorlnsl.%.4d",k);
		sv=fopen(fn,"rb");
		fread(xcorlnslt,sizeof(DP),(zd*xdv*3),sv);
		fclose(sv);
*/		
		sprintf(fn,"ycorl.%.4d",k);
		sv=fopen(fn,"rb");
		fread(ycorlt,sizeof(DP),(zd*yd*3),sv);
		fclose(sv);

		sprintf(fn,"uavg.%d",k);
		sv=fopen(fn,"rb");
		fread(uavgt,sizeof(DP),zd,sv);
		fclose(sv);

		for (a=0;a<(zd*3*xdv);a++)
		{
			*(xcorl+a) += (*(xcorlt+a));
			*(xcorlsl+a) += (*(xcorlslt+a));
			*(xcorlnsl+a) += (*(xcorlnslt+a));
		}
		for (a=0;a<(zd*3*yd);a++)
		{
			*(ycorl+a) += (*(ycorlt+a));
		}
		for (a=0;a<zd;a++)
			uavg[a] += uavgt[a];
	}

	for (a=0;a<(zd*3*xdv);a++)
	{
		*(xcorl+a) /= ((double)num);
		*(xcorlsl+a) /= ((double)num);
		*(xcorlnsl+a) /= ((double)num);
	}
	
	for (a=0;a<(zd*3*yd);a++)
	{
		*(ycorl+a) /= ((double)num);
	}

	for (k=0;k<zd;k++)
		uavg[k] /= ((double)num);
	
	for (k=0;k<up;k++)
	{
		for (j=0;j<yd;j++)
		{
			Fy(ycorl,j,k,0,yd,zd,3) = (Fy(ycorl,j,k,0,yd,zd,3)+Fy(ycorl,j,(zd-1-k),0,yd,zd,3))*0.5;
			Fy(ycorl,j,k,1,yd,zd,3) = (Fy(ycorl,j,k,1,yd,zd,3)+Fy(ycorl,j,(zd-1-k),1,yd,zd,3))*0.5;
			Fy(ycorl,j,k,2,yd,zd,3) = (Fy(ycorl,j,k,2,yd,zd,3)+Fy(ycorl,j,(zd-1-k),2,yd,zd,3))*0.5;
		}
	}
	
	for (k=0;k<up;k++)
	{
		for (i=0;i<xdv;i++)
		{
			Fx(xcorl,i,k,0,xdv,zd,3) = (Fx(xcorl,i,k,0,xdv,zd,3)+Fx(xcorl,i,(zd-1-k),0,xdv,zd,3))*0.5;
			Fx(xcorl,i,k,1,xdv,zd,3) = (Fx(xcorl,i,k,1,xdv,zd,3)+Fx(xcorl,i,(zd-1-k),1,xdv,zd,3))*0.5;
			Fx(xcorl,i,k,2,xdv,zd,3) = (Fx(xcorl,i,k,2,xdv,zd,3)+Fx(xcorl,i,(zd-1-k),2,xdv,zd,3))*0.5;

			Fx(xcorlsl,i,k,0,xdv,zd,3) = (Fx(xcorlsl,i,k,0,xdv,zd,3)+Fx(xcorlsl,i,(zd-1-k),0,xdv,zd,3))*0.5;
			Fx(xcorlsl,i,k,1,xdv,zd,3) = (Fx(xcorlsl,i,k,1,xdv,zd,3)+Fx(xcorlsl,i,(zd-1-k),1,xdv,zd,3))*0.5;
			Fx(xcorlsl,i,k,2,xdv,zd,3) = (Fx(xcorlsl,i,k,2,xdv,zd,3)+Fx(xcorlsl,i,(zd-1-k),2,xdv,zd,3))*0.5;

			Fx(xcorlnsl,i,k,0,xdv,zd,3) = (Fx(xcorlnsl,i,k,0,xdv,zd,3)+Fx(xcorlnsl,i,(zd-1-k),0,xdv,zd,3))*0.5;
			Fx(xcorlnsl,i,k,1,xdv,zd,3) = (Fx(xcorlnsl,i,k,1,xdv,zd,3)+Fx(xcorlnsl,i,(zd-1-k),1,xdv,zd,3))*0.5;
			Fx(xcorlnsl,i,k,2,xdv,zd,3) = (Fx(xcorlnsl,i,k,2,xdv,zd,3)+Fx(xcorlnsl,i,(zd-1-k),2,xdv,zd,3))*0.5;
		}
	}

	for (k=0;k<up;k++)
	{
		for (j=0;j<exty;j++)
		{
			Fy(ycorl,j,k,0,yd,zd,3) = (Fy(ycorl,j,k,0,yd,zd,3)+Fy(ycorl,(ydv-1-j),k,0,yd,zd,3))*0.5;
			Fy(ycorl,j,k,1,yd,zd,3) = (Fy(ycorl,j,k,1,yd,zd,3)+Fy(ycorl,(ydv-1-j),k,1,yd,zd,3))*0.5;
			Fy(ycorl,j,k,2,yd,zd,3) = (Fy(ycorl,j,k,2,yd,zd,3)+Fy(ycorl,(ydv-1-j),k,2,yd,zd,3))*0.5;
		}
	}
	
	for (k=0;k<up;k++)
	{
		for (i=0;i<extx;i++)
		{
			Fx(xcorl,i,k,0,xdv,zd,3) = (Fx(xcorl,i,k,0,xdv,zd,3)+Fx(xcorl,(xdv-1-i),k,0,xdv,zd,3))*0.5;
			Fx(xcorl,i,k,1,xdv,zd,3) = (Fx(xcorl,i,k,1,xdv,zd,3)+Fx(xcorl,(xdv-1-i),k,1,xdv,zd,3))*0.5;
			Fx(xcorl,i,k,2,xdv,zd,3) = (Fx(xcorl,i,k,2,xdv,zd,3)+Fx(xcorl,(xdv-1-i),k,2,xdv,zd,3))*0.5;

			Fx(xcorlsl,i,k,0,xdv,zd,3) = (Fx(xcorlsl,i,k,0,xdv,zd,3)+Fx(xcorlsl,(xdv-1-i),k,0,xdv,zd,3))*0.5;
			Fx(xcorlsl,i,k,1,xdv,zd,3) = (Fx(xcorlsl,i,k,1,xdv,zd,3)+Fx(xcorlsl,(xdv-1-i),k,1,xdv,zd,3))*0.5;
			Fx(xcorlsl,i,k,2,xdv,zd,3) = (Fx(xcorlsl,i,k,2,xdv,zd,3)+Fx(xcorlsl,(xdv-1-i),k,2,xdv,zd,3))*0.5;

			Fx(xcorlnsl,i,k,0,xdv,zd,3) = (Fx(xcorlnsl,i,k,0,xdv,zd,3)+Fx(xcorlnsl,(xdv-1-i),k,0,xdv,zd,3))*0.5;
			Fx(xcorlnsl,i,k,1,xdv,zd,3) = (Fx(xcorlnsl,i,k,1,xdv,zd,3)+Fx(xcorlnsl,(xdv-1-i),k,1,xdv,zd,3))*0.5;
			Fx(xcorlnsl,i,k,2,xdv,zd,3) = (Fx(xcorlnsl,i,k,2,xdv,zd,3)+Fx(xcorlnsl,(xdv-1-i),k,2,xdv,zd,3))*0.5;
		}
	}
 	
	sprintf(fn,"RX_z%d.dat",((int)ZPTGT));
#ifdef ZOH
	sprintf(fn,"RX_z%.2f.dat",(ZOH));
#endif
	sv=fopen(fn,"w");
	for (i=0;i<extx;i++)
	{
	  fprintf(sv,"%f ",(i+0.5)*utauavg/nu);
	  fprintf(sv,"%f ",Fx(xcorl,i,nz,0,xdv,zd,3)/Fx(xcorl,0,nz,0,xdv,zd,3));
	  fprintf(sv,"%f ",Fx(xcorl,i,nz,1,xdv,zd,3)/Fx(xcorl,0,nz,1,xdv,zd,3));
	  fprintf(sv,"%f\n",Fx(xcorl,i,nz,2,xdv,zd,3)/Fx(xcorl,0,nz,2,xdv,zd,3));
	}
	fclose(sv);
/*
	sprintf(fn,"RX_z%d-sl.dat",((int)ZPTGT));
#ifdef ZOH
	sprintf(fn,"RX_z%.2f-sl.dat",(ZOH));
#endif
	sv=fopen(fn,"w");
	for (i=0;i<extx;i++)
	{
	  fprintf(sv,"%f ",(i+0.5)*utauavg/nu);
	  fprintf(sv,"%f ",Fx(xcorlsl,i,nz,0,xdv,zd,3)/Fx(xcorlsl,0,nz,0,xdv,zd,3));
	  fprintf(sv,"%f ",Fx(xcorlsl,i,nz,1,xdv,zd,3)/Fx(xcorlsl,0,nz,1,xdv,zd,3));
	  fprintf(sv,"%f\n",Fx(xcorlsl,i,nz,2,xdv,zd,3)/Fx(xcorlsl,0,nz,2,xdv,zd,3));
	}
	fclose(sv);

	sprintf(fn,"RX_z%d-nsl.dat",((int)ZPTGT));
#ifdef ZOH
	sprintf(fn,"RX_z%.2f-nsl.dat",(ZOH));
#endif
	sv=fopen(fn,"w");
	for (i=0;i<extx;i++)
	{
	  fprintf(sv,"%f ",(i+0.5)*utauavg/nu);
	  fprintf(sv,"%f ",Fx(xcorlnsl,i,nz,0,xdv,zd,3)/Fx(xcorlnsl,0,nz,0,xdv,zd,3));
	  fprintf(sv,"%f ",Fx(xcorlnsl,i,nz,1,xdv,zd,3)/Fx(xcorlnsl,0,nz,1,xdv,zd,3));
	  fprintf(sv,"%f\n",Fx(xcorlnsl,i,nz,2,xdv,zd,3)/Fx(xcorlnsl,0,nz,2,xdv,zd,3));
	}
	fclose(sv);
*/	
	sprintf(fn,"RY_z%d.dat",((int)ZPTGT));
#ifdef ZOH
	sprintf(fn,"RY_z%.2f.dat",(ZOH));
#endif
	sv=fopen(fn,"w");
	for (j=0;j<exty;j++)
	{
	  fprintf(sv,"%f ",(j+0.5)*utauavg/nu);
	  fprintf(sv,"%f ",Fy(ycorl,j,nz,0,yd,zd,3)/Fy(ycorl,0,nz,0,yd,zd,3));
	  fprintf(sv,"%f ",Fy(ycorl,j,nz,1,yd,zd,3)/Fy(ycorl,0,nz,1,yd,zd,3));
	  fprintf(sv,"%f\n",Fy(ycorl,j,nz,2,yd,zd,3)/Fy(ycorl,0,nz,2,yd,zd,3));
	}
	fclose(sv);

	sprintf(fn,"uavg.dat");
	sv=fopen(fn,"w");
	for (k=0;k<up;k++)
		fprintf(sv,"%f %.9f\n",(*(z+k)+1)*Ret,0.5*(uavg[k]+uavg[zd-1-k])/utauavg);
	fclose(sv);
	
	free(ycorl);
	free(xcorl);
	free(xcorlsl);
	free(xcorlnsl);
	free(ycorlt);
	free(xcorlt);
	free(xcorlslt);
	free(xcorlnslt);
}	

int main()
{
	int xd,yd,zd,TM;
	int i,j,k,q,a,b,c,nut;
	int tagl=1,tagr=2;
	int cntr;
	int ss,up;
	FILE *epp;
	const DP C=1.;
	DP dx=1.;
	DP Fx=0.,Fy=0.,Fz=0.;
	DP Gx,Gy,Gz;
	DP tau=0.;
	DP rhozero=1.,rhozeroinv;
	DP t1,t2,t3,t4,t5;
	char ch;
	FILE *mysave,*tm;
	FILE *sv;
	char fn[20],fn2[20];
	DP utime,umax,dt;
	DP utavg;

	zd=zdv+2;
	xd=xdv;
	yd=ydv;

	correl_res(xd,yd,zd,TSTR,TEND);
	return(0);
}


	
	
	
	
	
	
	
	
	
	
