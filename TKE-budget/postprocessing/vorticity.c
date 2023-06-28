#include "definitions.h"

void Vorticity(int xd, int yd, int zd, int st, int en)
{
	int i,j,k,a;
	FILE *sv;
	FILE *tcp,*wstr;
	DP *mox,*moy,*moz,*fox,*foy,*foz,*moxt,*moyt,*mozt,*foxt,*foyt,*fozt;
	char fn[20];
	int num=(en-st)+1;
	DP *z,dplus,Ret;
	DP utauavg=0.,utvt,time,nu;
	int cnt=0,ct=0,rd;
	int up;
	DP ubavg=0.,ubt=0.,twall,val,tau;
	char stri[60];
	DP uttemp,uttemp1,ttemp;
	DP *du, *dudx;

	up=1+(zd-1)/2;
	num/=DT;
	z=(DP *)calloc(zd,sizeof(DP));

	fox = (DP *)calloc(zd,sizeof(DP));
	foy = (DP *)calloc(zd,sizeof(DP));
	foz = (DP *)calloc(zd,sizeof(DP));

	foxt = (DP *)calloc(zd,sizeof(DP));
	foyt = (DP *)calloc(zd,sizeof(DP));
	fozt = (DP *)calloc(zd,sizeof(DP));


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
	printf("In vorticity: ut=%f nu=%f\n",utauavg,nu);

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

	for (k=st;k<(en+1);k+=DT)
	{
		if (SLPF==0)
			sprintf(fn,"vort-fluc.%.3d.%.4d",0,k);
		else if (SLPF==-1)
			sprintf(fn,"vort-fluc-no-slip.%.3d.%.4d",0,k);
		else if (SLPF==-2)
			sprintf(fn,"vort-fluc-slip.%.3d.%.4d",0,k);
		else if (SLPF==-3)
			sprintf(fn,"vort-fluc-2d.%.3d.%.4d",0,k);
		sv=fopen(fn,"rb");
		fread(foxt,sizeof(DP),zd,sv);
		fread(foyt,sizeof(DP),zd,sv);
		fread(fozt,sizeof(DP),zd,sv);
		fclose(sv);

		for (a=0;a<zd;a++)
		{
			*(fox+a) += (*(foxt+a))*(*(foxt+a));
			*(foy+a) += (*(foyt+a))*(*(foyt+a));
			*(foz+a) += (*(fozt+a))*(*(fozt+a));
		}
	}

	for (k=0;k<zd;k++)
	{
		*(fox+k) /= num;
		*(foy+k) /= num;
		*(foz+k) /= num;
	}

	sprintf(stri,TEXT);
//	Ret=1.;

	sprintf(fn,"omega-rms-x-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(fox+k))+(*(fox+zd-1-k))))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-x-h-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(fox+k))+(*(fox+zd-1-k))))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
        fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(fox+k))+(*(fox+zd-1-k))))*nu/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"omega-rms-y-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(foy+k))+(*(foy+zd-1-k))))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-y-h-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(foy+k))+(*(foy+zd-1-k))))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
        fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(foy+k))+(*(foy+zd-1-k))))*nu/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"omega-rms-z-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(foz+k))+(*(foz+zd-1-k))))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-z-h-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(foz+k))+(*(foz+zd-1-k))))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
        fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(foz+k))+(*(foz+zd-1-k))))*nu/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	free(fox);
	free(foy);
	free(foz);

	free(foxt);
	free(foyt);
	free(fozt);
}