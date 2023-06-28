#include "definitions.h"

void KE_contour(int xd, int yd, int zd, int st, int en)
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
	double uu[yd][zd],vv[yd][zd],ww[yd][zd];
	double uut[yd][zd],vvt[yd][zd],wwt[yd][zd];
	double KE[yd][zd];

	up=1+(zd-1)/2;
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
	printf("Vorticity 2D: ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);

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

	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			uu[j][k] = 0.;
			vv[j][k] = 0.;
			ww[j][k] = 0.;

			uut[j][k] = 0.;
			vvt[j][k] = 0.;
			wwt[j][k] = 0.;
		}
	}

	for (k=st;k<(en+1);k+=DT)
	{
		if (SLPF==0)
                        sprintf(fn,"vel-field.%.3d.%.4d",0,k);
                else if (SLPF==-1)
                        sprintf(fn,"vel-field-no-slip.%.3d.%.4d",0,k);
                else if (SLPF==-2)
                        sprintf(fn,"vel-field-slip.%.3d.%.4d",0,k);
                else if (SLPF==-3)
                        sprintf(fn,"vel-field-2d.%.3d.%.4d",0,k);
                sv=fopen(fn,"rb");
                fread(&uut[0][0],sizeof(DP),(yd*zd),sv);
                fread(&vvt[0][0],sizeof(DP),(yd*zd),sv);
                fread(&wwt[0][0],sizeof(DP),(yd*zd),sv);
                fclose(sv);

		for (j=0;j<yd;j++)
		{
			for (a=1;a<(zd-1);a++)
			{
				uu[j][a] += uut[j][a];
				vv[j][a] += vvt[j][a];
				ww[j][a] += wwt[j][a];
			}
		}
	}

	for (j=0;j<yd;j++)
	{	
		for (k=0;k<zd;k++)
		{
			uu[j][k] /= ((double)num);
			vv[j][k] /= ((double)num);
			ww[j][k] /= ((double)num);
		}
	}

	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{	
			KE[j][k] = 0.5*(uu[j][k]+vv[j][k]+ww[j][k]);
		}
	}

	sprintf(stri,TEXT);
	sprintf(fn,"contour-k-plus-plus-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"ke contour\"\n");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\", \"z<sup>+</sup>\", \"k<sup>+</sup>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f\n",(j+0.5)*utauavg/nu,(k-0.5)*utauavg/nu,KE[j][k]/(utauavg*utauavg));
                }
        }
        fclose(tcp);

	sprintf(fn,"contour-k-bulk-plus-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"ke contour\"\n");
        fprintf(tcp,"VARIABLES = \"y/g\", \"z/g\", \"k<sup>+</sup>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f\n",(j+0.5)/gplus,(k-0.5)/gplus,KE[j][k]/(utauavg*utauavg));
                }
        }
        fclose(tcp);

	sprintf(fn,"contour-k-bulk-bulk-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"ke contour\"\n");
        fprintf(tcp,"VARIABLES = \"y/g\", \"z/g\", \"k<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f\n",(j+0.5)/gplus,(k-0.5)/gplus,KE[j][k]/(ub*ub));
                }
        }
        fclose(tcp);

	free(z);
}