#include "definitions.h"

void uvfl_pl(int xd, int yd, int zd, int tst)
{
	int i,j,k,a;
	FILE *sv;
	FILE *tcp,*wstr;
	char fn[20];
	DP *z,dplus,Ret;
	DP utauavg=0.,utvt,time,nu;
	int cnt=0,ct=0,rd;
	int up;
	DP ubavg=0.,ubt=0.,twall,val,tau,ub=0.;
	char stri[60];
	DP uttemp,uttemp1,ttemp;
	double ufl[xdv][ydv],vfl[xdv][ydv];

	up=1+(zd-1)/2;
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
	printf("uvflc: ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);
	
	sprintf(fn,"ufl4.%.4d",tst);
	sv=fopen(fn,"rb");
	fread(&ufl[0][0],xdv*ydv,sizeof(DP),sv);
	fclose(sv);

	sprintf(fn,"vfl4.%.4d",tst);
	sv=fopen(fn,"rb");
	fread(&vfl[0][0],xdv*ydv,sizeof(DP),sv);
	fclose(sv);

	sprintf(fn,"fluc.dat");
	sv=fopen(fn,"w");
	fprintf(sv,"TITLE = \"Cross flow contour\"\n");
	fprintf(sv,"VARIABLES = \"X<sup>+</sup>\", \"Y<sup>+</sup>\", \"u<sup>+</sup><sub>x</sub>\", \"u<sup>+</sup><sub>y</sub>\"\n");
	fprintf(sv,"ZONE I=%d, J=%d, F=POINT\n",xdv,ydv);

	for (j=0;j<ydv;j++)
	{
		for (i=0;i<xdv;i++)
		{
			fprintf(sv,"%f %f %f %f\n",(i+0.5)*utauavg/nu,(j+0.5)*utauavg/nu,ufl[i][j]/utauavg,vfl[i][j]/utauavg);
		}
	}
	fclose(sv);

}