#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* postporcessing the results of energy spectra code */
#define xdv 512           //Number of grid points in x
#define ydv 256           //Number of grid points in y
#define zdv 221           //Number of grid points in z

#define Reb 7200.0000     //bulk Reynolds number
#define ubulk 0.0866034   //nominal bulk velocity
#define TSTR 0            //begining time
#define TEND 1            //end time
#define DT 1              //time step for calculating the averages
#define ZPTGT 10          //approximate distance from the walls in wall units

#define ALPHA 1.35        //physical channel size in x, Lx=2xPi/alpha
#define BETA 2.5          //physical channel size in y, Ly=2xPi/beta

typedef double DP;
#define PI (4.*atan(1.))

void e_spact(int xd, int yd, int zd, int st, int en)
{
	int i,j,k,a,b;
	double Euux[xd][zd],Evvx[xd][zd],Ewwx[xd][zd];
	double Euuxt[xd][zd],Evvxt[xd][zd],Ewwxt[xd][zd];
	double Euuxe[xd/2],Evvxe[xd/2],Ewwxe[xd/2];
	double urms,urms2,vrms,vrms2,wrms,wrms2;
	
	double Euuy[yd/2 +1][zd],Evvy[yd/2 +1][zd],Ewwy[yd/2 +1][zd];
	double Euuyt[yd/2 +1][zd],Evvyt[yd/2 +1][zd],Ewwyt[yd/2 +1][zd];
	double Euuye[yd/2],Evvye[yd/2],Ewwye[yd/2];
	
	double Kx[xd/2],Ky[yd/2];
	double z[zd],zp,zn;
	int nz;
	
	double fac;
	
	int num=(en-st)+1;
	int up=1+(zd-1)/2;
	
	double nu,tau,utauavg,uttemp,uttemp1,ttemp,dplus,Ret;
	int ct=0;
	FILE *wstr,*fh;
	char fn[80];
	int nxs,nys;
	
	up=1+(zd-1)/2;
	num /= DT;
	
	for (i=0;i<xd/2;i++)
	  Kx[i] = i*ALPHA/(0.5*(zd-2.));
	
	for (j=0;j<yd/2;j++)
	  Ky[j] = j*BETA/(0.5*(zd-2.));
	
	//p-grad.txt contains the time history of pressure gradient during the simulation
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
	printf("ut=%f nu=%f  Ret=%f\n",utauavg,nu,utauavg*0.5*(zd-2)/nu);
	
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
	
	/* Initialization */
	printf("initialization\n");
	for (a=0;a<xd;a++)
	{
		for (b=0;b<zd;b++)
		{
			Euux[a][b] = 0.;
			Evvx[a][b] = 0.;
			Ewwx[a][b] = 0.;
		}
	}
	for (a=0;a<(yd/2 +1);a++)
	{
		for (b=0;b<zd;b++)
		{
			Euuy[a][b] = 0.;
			Evvy[a][b] = 0.;
			Ewwy[a][b] = 0.;
		}
	}
	for (a=0;a<xd/2;a++)
	{
		Euuxe[a] = 0.;
		Evvxe[a] = 0.;
		Ewwxe[a] = 0.;
	}
	for (a=0;a<yd/2;a++)
	{
		Euuye[a] = 0.;
		Evvye[a] = 0.;
		Ewwye[a] = 0.;
	}
	
	printf("reading\n");
	/* nz is the z plane of interest */
	for (k=st;k<(en+1);k+=DT)
	{
		sprintf(fn,"Euux.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Euuxt[0][0],sizeof(DP),(zd*xd),fh);
		fclose(fh);
		
		sprintf(fn,"Evvx.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Evvxt[0][0],sizeof(DP),(zd*xd),fh);
		fclose(fh);
		
		sprintf(fn,"Ewwx.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Ewwxt[0][0],sizeof(DP),(zd*xd),fh);
		fclose(fh);
		
		sprintf(fn,"Euuy.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Euuyt[0][0],sizeof(DP),(zd*(yd/2 +1)),fh);
		fclose(fh);
		
		sprintf(fn,"Evvy.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Evvyt[0][0],sizeof(DP),(zd*(yd/2 +1)),fh);
		fclose(fh);
		
		sprintf(fn,"Ewwy.%.4d",k);
		fh=fopen(fn,"rb");
		fread(&Ewwyt[0][0],sizeof(DP),(zd*(yd/2 +1)),fh);
		fclose(fh);
//		printf("adding\n");
		for (a=0;a<xd;a++)
		{
			for (b=0;b<zd;b++)
			{
				Euux[a][b] += Euuxt[a][b];
				Evvx[a][b] += Evvxt[a][b];
				Ewwx[a][b] += Ewwxt[a][b];
			}
		}
		
		for (a=0;a<(yd/2 +1);a++)
		{
			for (b=0;b<zd;b++)
			{
				Euuy[a][b] += Euuyt[a][b];
				Evvy[a][b] += Evvyt[a][b];
				Ewwy[a][b] += Ewwyt[a][b];
			}
		}
	}
	printf("end of reading \n");
	for (a=0;a<xd;a++)
	{
		for (b=0;b<zd;b++)
		{
			Euux[a][b] /= num;
			Evvx[a][b] /= num;
			Ewwx[a][b] /= num;
		}
	}
	for (a=0;a<(yd/2 +1);a++)
	{
		for (b=0;b<zd;b++)
		{
			Euuy[a][b] /= num;
			Evvy[a][b] /= num;
			Ewwy[a][b] /= num;
		}
	}

	urms=urms2=0.;
        vrms=vrms2=0.;
        wrms=wrms2=0.;
        for (a=0;a<xd;a++)
        {
		urms += 0.5*(Euux[a][nz]+Euux[a][zd-1-nz]);
                vrms += 0.5*(Evvx[a][nz]+Evvx[a][zd-1-nz]);
                wrms += 0.5*(Ewwx[a][nz]+Ewwx[a][zd-1-nz]);
	}
        for (a=0;a<(yd/2 +1);a++)
        {
                urms2 += 0.5*(Euuy[a][nz]+Euuy[a][zd-1-nz]);
                vrms2 += 0.5*(Evvy[a][nz]+Evvy[a][zd-1-nz]);
                wrms2 += 0.5*(Ewwy[a][nz]+Ewwy[a][zd-1-nz]);
        }

	printf("urms =%.12f,%f vrms =%.12f,%f wrms =%.12f,%f\n",urms,sqrt(urms)/utauavg,vrms,sqrt(vrms)/utauavg,wrms,sqrt(wrms)/utauavg);
        printf("urms2=%.12f,%f vrms2=%.12f,%f wrms2=%.12f,%f\n",urms2,sqrt(urms2)/utauavg,vrms2,sqrt(vrms2)/utauavg,wrms2,sqrt(wrms2)/utauavg);
	
	/* X spectra */
	a=0;
	Euuxe[a] = 0.5*(Euux[a][nz] + Euux[a][zd-1-nz]);
	Evvxe[a] = 0.5*(Evvx[a][nz] + Evvx[a][zd-1-nz]);
	Ewwxe[a] = 0.5*(Ewwx[a][nz] + Ewwx[a][zd-1-nz]);

	for (a=1;a<xd/2;a++)
	{
		Euuxe[a] = 0.5*(Euux[a][nz]+Euux[a][zd-1-nz]);// + 0.5*(Euux[xd-a][nz]+Euux[xd-a][zd-1-nz]);
		Evvxe[a] = 0.5*(Evvx[a][nz]+Evvx[a][zd-1-nz]);// + 0.5*(Evvx[xd-a][nz]+Evvx[xd-a][zd-1-nz]);
		Ewwxe[a] = 0.5*(Ewwx[a][nz]+Ewwx[a][zd-1-nz]);// + 0.5*(Ewwx[xd-a][nz]+Ewwx[xd-a][zd-1-nz]);
	}

	
	/* Y spectra */
	for (a=0;a<yd/2;a++)
	{
		Euuye[a] = 0.5*(Euuy[a][nz]+Euuy[a][zd-1-nz]);
		Evvye[a] = 0.5*(Evvy[a][nz]+Evvy[a][zd-1-nz]);
		Ewwye[a] = 0.5*(Ewwy[a][nz]+Ewwy[a][zd-1-nz]);
	}

	/* Dealiasing */
	nxs = 2*(xd/2)/3;
	nys = 2*((yd/2)/3);
/*	for (a=1;a<xd/2;a++)
	{
		if (a > nxs)
		{
			Euuxe[a] = 0.;
			Evvxe[a] = 0.;
			Ewwxe[a] = 0.;
		}
	}
	for (a=0;a<yd/2;a++)
	{
		if (a > nys)
		{
			Euuye[a] = 0.;
			Evvye[a] = 0.;
			Ewwye[a] = 0.;
		}
	}
*/	
	sprintf(fn,"EX-zpl%d.dat",ZPTGT);
	fh=fopen(fn,"w");
	fprintf(fh,"# Kxplus Euuxplus Evvxpl Ewwxpl\n");
	for (a=0;a<xd/2;a++)
		fprintf(fh,"%.14f %.14f %.14f %.14f\n",Kx[a]*nu/utauavg,Euuxe[a]*0.5*(zd-2.)/(ALPHA*nu*utauavg),Evvxe[a]*0.5*(zd-2.)/(ALPHA*nu*utauavg),Ewwxe[a]*0.5*(zd-2.)/(ALPHA*nu*utauavg));
	fclose(fh);
	
	sprintf(fn,"EY-zpl%d.dat",ZPTGT);
	fh=fopen(fn,"w");
	fprintf(fh,"# Kyplus Euuyplus Evvypl Ewwypl\n");
	for (a=0;a<yd/2;a++)
		fprintf(fh,"%.14f %.14f %.14f %.14f\n",Ky[a]*nu/utauavg,Euuye[a]*0.5*(zd-2.)/(BETA*nu*utauavg),Evvye[a]*0.5*(zd-2.)/(BETA*nu*utauavg),Ewwye[a]*0.5*(zd-2.)/(BETA*nu*utauavg));
	fclose(fh);

	urms=urms2=0.;
	vrms=vrms2=0.;
	wrms=wrms2=0.;
	for (a=0;a<xd/2;a++)
	{
		urms += Euuxe[a];
		vrms += Evvxe[a];
		wrms += Ewwxe[a];
	}
	for (a=0;a<yd/2;a++)
        {
                urms2 += Euuye[a];
                vrms2 += Evvye[a];
                wrms2 += Ewwye[a];
        }
	printf("urms =%.12f,%f vrms =%.12f,%f wrms =%.12f,%f\n",urms,sqrt(urms)/utauavg,vrms,sqrt(vrms)/utauavg,wrms,sqrt(wrms)/utauavg);
	printf("urms2=%.12f,%f vrms2=%.12f,%f wrms2=%.12f,%f\n",urms2,sqrt(urms2)/utauavg,vrms2,sqrt(vrms2)/utauavg,wrms2,sqrt(wrms2)/utauavg);

}

int main()
{
	int i,j,k;
	int xd,yd,zd;
	
	xd=xdv;
	yd=ydv;
	zd=zdv+2;
	
	e_spact(xd, yd, zd, TSTR, TEND);
  
  return 0;
}
