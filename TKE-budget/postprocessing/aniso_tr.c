#include "definitions.h"

void aniso_tr(int xd, int yd, int zd, int st, int en)
{
	int i,j,k,a;
	FILE *sv;
	FILE *tcp,*wstr;
	char fn[20];
	int num=(en-st)+1;
	DP *z,dplus,Ret;
	DP utauavg=0.,utvt,time,nu;
	int cnt=0,ct=0,rd;
	int up;
	DP ubavg=0.,ubt=0.,twall,val,tau,ub=0.;
	char stri[60];
	DP uttemp,uttemp1,ttemp;
	DP *utin,*vtin,*wtin,*uw,*dtin;
	DP *utint,*vtint,*wtint,*uwt,*dtint;
	DP *uv,*uvt,*vw,*vwt;
	DP *UW,*UWt,*KE;
	DP *TMP;

	DP A[3][3]={0.};
	DP IIa[zd],IIIa[zd];

	utin = (DP *)calloc(zd,sizeof(DP));
	vtin = (DP *)calloc(zd,sizeof(DP));
	wtin = (DP *)calloc(zd,sizeof(DP));
	uw = (DP *)calloc(zd,sizeof(DP));

	uv = (DP *)calloc(zd,sizeof(DP));
	uvt = (DP *)calloc(zd,sizeof(DP));
	vw = (DP *)calloc(zd,sizeof(DP));
	vwt = (DP *)calloc(zd,sizeof(DP));

	UW = (DP *)calloc(zd,sizeof(DP));
	UWt= (DP *)calloc(zd,sizeof(DP));

	utint = (DP *)calloc(zd,sizeof(DP));
	vtint = (DP *)calloc(zd,sizeof(DP));
	wtint = (DP *)calloc(zd,sizeof(DP));
	uwt = (DP *)calloc(zd,sizeof(DP));

	TMP = (DP *)calloc(zd,sizeof(DP));
	KE = (DP *)calloc(zd,sizeof(DP));

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
	printf("aniso_tr: ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);

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
	
	for (k=st;k<(en+1);k+=DT)
	{
		if (SLPF==0)
			sprintf(fn,"turb-fld.%.3d.%.4d",0,k);
		else if (SLPF==-1)
			sprintf(fn,"turb-fld-no-slip.%.3d.%.4d",0,k);
		else if (SLPF==-2)
			sprintf(fn,"turb-fld-slip.%.3d.%.4d",0,k);
		else if (SLPF==-3)
			sprintf(fn,"turb-fld-2d.%.3d.%.4d",0,k);
		sv=fopen(fn,"rb");
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(utint,sizeof(DP),zd,sv);
		fread(vtint,sizeof(DP),zd,sv);
		fread(wtint,sizeof(DP),zd,sv);
		fread(uvt,sizeof(DP),zd,sv);
		fread(uwt,sizeof(DP),zd,sv);
		fread(vwt,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);
		fread(TMP,sizeof(DP),zd,sv);

		fread(UWt,sizeof(DP),zd,sv);
		fclose(sv);

		for (a=0;a<zd;a++)
		{
			*(uv+a) += *(uvt+a);
			*(vw+a) += *(vwt+a);

			*(utin+a) += (*(utint+a))*(*(utint+a));
			*(vtin+a) += (*(vtint+a))*(*(vtint+a));
			*(wtin+a) += (*(wtint+a))*(*(wtint+a));
			
			*(uw+a) += (*(uwt+a));
			*(UW+a) += (*(UWt+a));
		}
	}

	for (k=0;k<zd;k++)
	{
		*(uv+k) /= num;
		*(vw+k) /= num;
		*(UW+k) /= num;

		*(utin+k) /= num;
		*(vtin+k) /= num;
		*(wtin+k) /= num;
		*(uw+k) /= num;
	}

	for (k=0;k<zd;k++)
	{
		*(KE+k) = ((*(utin+k))+(*(vtin+k))+(*(wtin+k)));
		IIa[k] = 0.;
		IIIa[k] = 0.;
	}

	for (k=1;k<(zd-1);k++)
	{
		A[0][0] = (*(utin+k))/(*(KE+k))-1./3.;
		A[0][1] = (*(uv+k))/(*(KE+k));
		A[1][0] = A[0][1];
		A[0][2] = (*(uw+k))/(*(KE+k));
		A[2][0] = A[0][2];
		A[1][1] = (*(vtin+k))/(*(KE+k))-1./3.;
		A[1][2] = (*(vw+k))/(*(KE+k));
		A[2][1] = A[1][2];
		A[2][2] = (*(wtin+k))/(*(KE+k))-1./3.;

		for (i=0;i<3;i++)
		{
			for (j=0;j<3;j++)
			{
				IIa[k] += A[i][j]*A[j][i];
				for (a=0;a<3;a++)
				{
					IIIa[k] += A[i][j]*A[j][a]*A[a][i];
				}
			}
		}
	}

	sprintf(fn,"Aniso-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"II\",\"III\"\n");
	for (k=up-1;k>0;k--)
		fprintf(tcp,"%f %.9f %.9f\n",((*(z+k))+1)*Ret,IIa[k],IIIa[k]);
	fprintf(tcp,"\n");
	fclose(tcp);



	free(utin);
	free(vtin);
	free(wtin);

	free(uv);
	free(uw);
	free(vw);
	free(UW);

	free(KE);

	free(utint);
	free(vtint);
	free(wtint);

	free(uvt);
	free(uwt);
	free(vwt);
	free(UWt);

	free(TMP);
	free(z);
}