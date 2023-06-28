#include "definitions.h"

void budget_Line_results_collected(int xd, int yd, int zd, int st, int en)
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
	DP Pii[yd][zd],Epii[yd][zd],Trii[yd][zd],VdKdy[yd][zd],WdKdz[yd][zd];
	DP Piit[yd][zd],Epiit[yd][zd],Triit[yd][zd],VdKdyt[yd][zd],WdKdzt[yd][zd];
	DP P1[yd][zd],P2[yd][zd],P3[yd][zd],P4[yd][zd],P5[yd][zd],P6[yd][zd];
	DP P1t[yd][zd],P2t[yd][zd],P3t[yd][zd],P4t[yd][zd],P5t[yd][zd],P6t[yd][zd];
	DP Temp[yd][zd];
	DP U[yd][zd],Ut[yd][zd];
	DP uAVGt[yd][zd],uAVG[yd][zd];
	DP Uns[zd],Us[zd];
	DP UW[yd][zd],UWt[yd][zd],uw[yd][zd],uwt[yd][zd],uw1[zd],uw3[zd],UW3[zd],UW1[zd];
	DP Tiir[yd][zd],Tiis[yd][zd],Tiipai[yd][zd];

	DP Pii1[zd],Epii1[zd],Trii1[zd],VdKdy1[zd],WdKdz1[zd],U1[zd],Tiir1[zd],Tiipai1[zd],Tiis1[zd];
	DP Pii3[zd],Epii3[zd],Trii3[zd],VdKdy3[zd],WdKdz3[zd],U3[zd],Tiir3[zd],Tiipai3[zd],Tiis3[zd];
	
	DP P11[zd],P21[zd],P31[zd],P41[zd],P51[zd],P61[zd];
	DP P13[zd],P23[zd],P33[zd],P43[zd],P53[zd],P63[zd];

	DP S[yd],Prd[zd],Eps[zd];
	DP wstravg[yd],u1,u2,u3,Pii1np[zd],Epii1np[zd],Trii1np[zd],P11np[zd],P21np[zd],P31np[zd];

	double utaut_nsl[zd],utaut_sl[zd],utaut_avg[zd];
        DP uttotal,dsl[zd],dnsl[zd],dtot[zd];
	DP dudxsl[zd],dudxnsl[zd];
	
	for (j=0;j<yd;j++)
	{
		S[j]= 0.;
		wstravg[j] = 0.;
	}

	for (j=wst;j<yd;j+=(gplus+wplus))
	{
		for (a=0;a<gplus;a++)
			S[(j+a)]=-1;
	}

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
	printf("budget_Line_results_collected: ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);

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
	{	
		Pii1[k] = 0.;
		Epii1[k] = 0.;
		Trii1[k] = 0.;
		VdKdy1[k] = 0.;
		WdKdz1[k] = 0.;
		U1[k] = 0.;
		P11[k] = 0.;
		P21[k] = 0.;
		P31[k] = 0.;
		P41[k] = 0.;
		P51[k] = 0.;
		P61[k] = 0.;

		Pii3[k] = 0.;
		Epii3[k] = 0.;
		Trii3[k] = 0.;
		VdKdy3[k] = 0.;
		WdKdz3[k] = 0.;
		U3[k] = 0.;
		P13[k] = 0.;
		P23[k] = 0.;
		P33[k] = 0.;
		P43[k] = 0.;
		P53[k] = 0.;
		P63[k] = 0.;

		Prd[k]=0.;
		Eps[k]=0.;

		Uns[k] = 0.;
		Us[k] = 0.;

		uw1[k] = 0.;
		uw3[k] = 0.;

		UW1[k] = 0.;
		UW3[k] = 0.;

		Tiis1[k] = 0.;
		Tiir1[k] = 0.;
		Tiipai1[k] = 0.;

		Tiis3[k] = 0.;
                Tiir3[k] = 0.;
                Tiipai3[k] = 0.;
	}

	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			Pii[j][k] = 0.;
			Epii[j][k] = 0.;
			VdKdy[j][k] = 0.;
			WdKdz[j][k] = 0.;
			Trii[j][k] = 0.;

			P1[j][k] = 0.;
			P2[j][k] = 0.;
			P3[j][k] = 0.;
			P4[j][k] = 0.;
			P5[j][k] = 0.;
			P6[j][k] = 0.;

			uAVGt[j][k] = 0.;
			uAVG[j][k] = 0.;

			UW[j][k] = 0.;
			UWt[j][k] = 0.;
			uw[j][k] = 0.;
			uwt[j][k] = 0.;

			Tiir[j][k] = 0.;
			Tiipai[j][k] = 0.;
			Tiis[j][k] = 0.;
		}
	}

	for (k=0;k<zd;k++)
        {
                Trii1np[k] = 0.;
                Epii1np[k] = 0.;
                Pii1np[k] = 0.;

                P11np[k] = 0.;
                P21np[k] = 0.;
                P31np[k] = 0.;
        }
	
	for (k=st;k<(en+1);k+=DT)
	{
		if (SLPF==0)
			sprintf(fn,"Budget-2d-field.%.3d.%.4d",0,k);
		else if (SLPF==-1)
			sprintf(fn,"Budget-2d-field.%.3d.%.4d",0,k);
		else if (SLPF==-2)
			sprintf(fn,"Budget-2d-field.%.3d.%.4d",0,k);
		else if (SLPF==-3)
			sprintf(fn,"Budget-2d-field.%.3d.%.4d",0,k);
		sv=fopen(fn,"rb");
		fread(&Piit[0][0],sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(&Triit[0][0],sizeof(DP),(yd*zd),sv);
		for (j=0;j<yd;j++)
		{
			for (a=0;a<zd;a++)
			{
				Trii[j][a] += Triit[j][a];
				Tiir[j][a] += Triit[j][a];
			}
		}
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(&Triit[0][0],sizeof(DP),(yd*zd),sv);
		for (j=0;j<yd;j++)
		{
			for (a=0;a<zd;a++)
			{
				Trii[j][a] += Triit[j][a];
				Tiipai[j][a] += Triit[j][a];
			}
		}
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(&Triit[0][0],sizeof(DP),(yd*zd),sv);
		for (j=0;j<yd;j++)
		{
			for (a=0;a<zd;a++)
			{
				Trii[j][a] += Triit[j][a];
				Tiis[j][a] += Triit[j][a];
			}
		}
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(&Epiit[0][0],sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(&VdKdyt[0][0],sizeof(DP),(yd*zd),sv);
		fread(&WdKdzt[0][0],sizeof(DP),(yd*zd),sv);
		fclose(sv);

		if (SLPF==0)
			sprintf(fn,"TP-Breakdown-2d-field.%.3d.%.4d",0,k);
		else if (SLPF==-1)
			sprintf(fn,"TP-Breakdown-2d-field.%.3d.%.4d",0,k);
		else if (SLPF==-2)
			sprintf(fn,"TP-Breakdown-2d-field.%.3d.%.4d",0,k);
		else if (SLPF==-3)
			sprintf(fn,"TP-Breakdown-2d-field.%.3d.%.4d",0,k);
		sv=fopen(fn,"rb");
		fread(&P1t[0][0],sizeof(DP),(yd*zd),sv);
		fread(&P2t[0][0],sizeof(DP),(yd*zd),sv);
		fread(&P3t[0][0],sizeof(DP),(yd*zd),sv);
		fread(&P4t[0][0],sizeof(DP),(yd*zd),sv);
		fread(&P5t[0][0],sizeof(DP),(yd*zd),sv);
		fread(&P6t[0][0],sizeof(DP),(yd*zd),sv);
		fclose(sv);

		sprintf(fn,"uAVG-2d.%.3d.%.4d",0,k);
                sv=fopen(fn,"rb");
                fread(&uAVGt[0][0],sizeof(DP),yd*zd,sv);
                fclose(sv);
		
		sprintf(fn,"vel-field-2d.%.3d.%.4d",0,k);
                sv=fopen(fn,"rb");
		fseek(sv,4*yd*zd*sizeof(double),SEEK_SET);
                fread(&uwt[0][0],sizeof(DP),(yd*zd),sv);
		fseek(sv,6*yd*zd*sizeof(double),SEEK_SET);
		fread(&UWt[0][0],sizeof(DP),(zd*yd),sv);
                fclose(sv);

		for (j=0;j<yd;j++)
		{
			for (a=1;a<(zd-1);a++)
			{
				Pii[j][a] += Piit[j][a];
				Epii[j][a] += Epiit[j][a];
				VdKdy[j][a] += VdKdyt[j][a];
				WdKdz[j][a] += WdKdzt[j][a];

				P1[j][a] += P1t[j][a];
				P2[j][a] += P2t[j][a];
				P3[j][a] += P3t[j][a];
				P4[j][a] += P4t[j][a];
				P5[j][a] += P5t[j][a];
				P6[j][a] += P6t[j][a];

				uAVG[j][a] += uAVGt[j][a];

				uw[j][a] += uwt[j][a];
                                UW[j][a] += UWt[j][a];
			}
		}
	}

	for (j=0;j<yd;j++)
	{	
		for (k=0;k<zd;k++)
		{
			Trii[j][k] /= ((double)num);
			Pii[j][k] /= ((double)num);
			Epii[j][k] /= ((double)num);
			VdKdy[j][k] /= ((double)num);
			WdKdz[j][k] /= ((double)num);

			P1[j][k] /= ((double)num);
			P2[j][k] /= ((double)num);
			P3[j][k] /= ((double)num);
			P4[j][k] /= ((double)num);
			P5[j][k] /= ((double)num);
			P6[j][k] /= ((double)num);

			uAVG[j][k] /= ((double)num);

			uw[j][k] /= ((double)num);
                        UW[j][k] /= ((double)num);
	
			Tiir[j][k] /= ((double)num);
			Tiipai[j][k] /= ((double)num);
			Tiis[j][k] /= ((double)num);
		}
	}
	ct=0;
	for (j=0;j<yd;j++)
	{
		if (S[j]==-1)
		{
			for (k=0;k<zd;k++)
			{
				Trii3[k] += Trii[j][k];
				Epii3[k] += Epii[j][k];
				Pii3[k] += Pii[j][k];
				VdKdy3[k] += VdKdy[j][k];
				WdKdz3[k] += WdKdz[j][k];
				P13[k] += P1[j][k];
				P23[k] += P2[j][k];
				P33[k] += P3[j][k];
				P43[k] += P4[j][k];
				P53[k] += P5[j][k];
				P63[k] += P6[j][k];
		
				Us[k] += uAVG[j][k];
				uw3[k] += uw[j][k];
				UW3[k] += UW[j][k];

				Tiir3[k] += Tiir[j][k];
				Tiipai3[k] += Tiipai[j][k];
				Tiis3[k] += Tiis[j][k];
			}
			ct++;
		}
	}
	printf("c = %d\n",ct);
/*	ct=0;
	for (j=wst;j<yd;j+=(gplus+wplus))
	{
		for (a=0;a<gplus;a++)
		{
			for (k=0;k<zd;k++)
			{
				Trii3[k] += Trii[j+a][k];
				Epii3[k] += Epii[j+a][k];
				Pii3[k] += Pii[j+a][k];
				VdKdy3[k] += VdKdy[j+a][k];
				WdKdz3[k] += WdKdz[j+a][k];
				P13[k] += P1[j+a][k];
				P23[k] += P2[j+a][k];
				P33[k] += P3[j+a][k];
				P43[k] += P4[j+a][k];
				P53[k] += P5[j+a][k];
				P63[k] += P6[j+a][k];
			}
			ct++;
		}
	}
*/
	for (k=0;k<zd;k++)
	{
		Trii3[k] /= ((double)ct);
		Epii3[k] /= ((double)ct);
		Pii3[k] /= ((double)ct);
		VdKdy3[k] /= ((double)ct);
		WdKdz3[k] /= ((double)ct);
		P13[k] /= ((double)ct);
		P23[k] /= ((double)ct);
		P33[k] /= ((double)ct);
		P43[k] /= ((double)ct);
		P53[k] /= ((double)ct);
		P63[k] /= ((double)ct);

		Us[k] /= ((double)ct);
		uw3[k] /= ((double)ct);
		UW3[k] /= ((double)ct);

		Tiir3[k] /= ((double)ct);
                Tiipai3[k] /= ((double)ct);
                Tiis3[k] /= ((double)ct);
	}

	ct=0;
	for (j=0;j<yd;j++)
	{
		if (S[j]==0)
		{
			for (k=0;k<zd;k++)
			{
				Trii1[k] += Trii[j][k];
				Epii1[k] += Epii[j][k];
				Pii1[k] += Pii[j][k];
				VdKdy1[k] += VdKdy[j][k];
				WdKdz1[k] += WdKdz[j][k];
				P11[k] += P1[j][k];
				P21[k] += P2[j][k];
				P31[k] += P3[j][k];
				P41[k] += P4[j][k];
				P51[k] += P5[j][k];
				P61[k] += P6[j][k];

				Uns[k] += uAVG[j][k];
				uw1[k] += uw[j][k];
				UW1[k] += UW[j][k];

				Tiir1[k] += Tiir[j][k];
				Tiipai1[k] += Tiipai[j][k];
				Tiis1[k] += Tiis[j][k];
			}
			ct++;
		}
	}
	printf("c = %d\n",ct);
/*
	ct=0;
	for (j=(wst+gplus);j<(yd-wst);j+=(gplus+wplus))
	{
		for (a=0;a<wplus;a++)
		{
			for (k=0;k<zd;k++)
			{
				Trii1[k] += Trii[j+a][k];
				Epii1[k] += Epii[j+a][k];
				Pii1[k] += Pii[j+a][k];
				VdKdy1[k] += VdKdy[j+a][k];
				WdKdz1[k] += WdKdz[j+a][k];
				P11[k] += P1[j+a][k];
				P21[k] += P2[j+a][k];
				P31[k] += P3[j+a][k];
				P41[k] += P4[j+a][k];
				P51[k] += P5[j+a][k];
				P61[k] += P6[j+a][k];
			}
			ct++;
		}
	}
	for (a=0;a<wst;a++)
	{
		for (k=0;k<zd;k++)
		{
			Trii1[k] += Trii[a][k];
			Epii1[k] += Epii[a][k];
			Pii1[k] += Pii[a][k];
			VdKdy1[k] += VdKdy[a][k];
			WdKdz1[k] += WdKdz[a][k];
			P11[k] += P1[a][k];
			P21[k] += P2[a][k];
			P31[k] += P3[a][k];
			P41[k] += P4[a][k];
			P51[k] += P5[a][k];
			P61[k] += P6[a][k];

			Trii1[k] += Trii[yd-wst+a][k];
			Epii1[k] += Epii[yd-wst+a][k];
			Pii1[k] += Pii[yd-wst+a][k];
			VdKdy1[k] += VdKdy[yd-wst+a][k];
			WdKdz1[k] += WdKdz[yd-wst+a][k];
			P11[k] += P1[yd-wst+a][k];
			P21[k] += P2[yd-wst+a][k];
			P31[k] += P3[yd-wst+a][k];
			P41[k] += P4[yd-wst+a][k];
			P51[k] += P5[yd-wst+a][k];
			P61[k] += P6[yd-wst+a][k];
		}
		ct++;
	}
*/
	for (k=0;k<zd;k++)
	{
		Trii1[k] /= ((double)ct);
		Epii1[k] /= ((double)ct);
		Pii1[k] /= ((double)ct);
		VdKdy1[k] /= ((double)ct);
		WdKdz1[k] /= ((double)ct);
		P11[k] /= ((double)ct);
		P21[k] /= ((double)ct);
		P31[k] /= ((double)ct);
		P41[k] /= ((double)ct);
		P51[k] /= ((double)ct);
		P61[k] /= ((double)ct);

		Uns[k] /= ((double)ct);

		uw1[k] /= ((double)ct);
		UW1[k] /= ((double)ct);

		Tiir1[k] /= ((double)ct);
		Tiipai1[k] /= ((double)ct);
		Tiis1[k] /= ((double)ct);
	}

	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			Eps[k] += Epii[j][k];
			Prd[k] += Pii[j][k];
		}
	}
	for (k=0;k<zd;k++)
	{
		Prd[k] /= ((double)yd);
		Eps[k] /= ((double)yd);
	}
//****************************************************************************************
	for (j=0;j<yd;j++)
                wstravg[j]=0.;

        for (j=0;j<yd;j++)
        {
                  u1=uAVG[j][1];
                  u2=uAVG[j][2];
                  u3=uAVG[j][3];

                  wstravg[j]=nu*(-2.*u1+3.*u2-u3);
        }

        for (j=wst;j<yd;j+=(gplus+wplus))
        {
                for (a=0;a<gplus;a++)
                {
                        wstravg[j+a] = 0.;
                }
        }

        for (j=wst;j<yd;j+=(gplus+wplus))
        {
                for (a=0;a<gplus;a++)
                {
                        wstravg[j+a] = 0.;
                }
        }

	for (k=3;k<up;k++)
        {
                dudxsl[k] = (2./3.)*(Us[k+1]-Us[k-1])-(1./12.)*(Us[k+2]-Us[k-2]);
                dudxnsl[k] = (2./3.)*(Uns[k+1]-Uns[k-1])-(1./12.)*(Uns[k+2]-Uns[k-2]);
        }

        dudxsl[0] = 0.;
        dudxnsl[0] = 0.;

        dudxnsl[0] = 2.*Uns[1];

	dudxsl[1] = 0.5*(Us[2]-Us[1]);
        dudxnsl[1] = 0.5*(Uns[2]+Uns[1]);

        dudxsl[2] = 0.5*(Us[3]-Us[1]);
        dudxnsl[2] = 0.5*(Uns[3]-Uns[1]);

	for (k=0;k<(up-1);k++)
        {
                utaut_nsl[k] = sqrt((nu*dudxnsl[k]-0.5*(uw1[k]-uw1[zd-1-k]))/(1.-((*(z+k))+1)));
                utaut_sl[k] = sqrt((nu*dudxsl[k] - 0.5*(uw3[k]-uw3[zd-1-k]))/(1.-((*(z+k))+1)));
                utaut_avg[k] = sqrt((gplus*utaut_sl[k]*utaut_sl[k]+wplus*utaut_nsl[k]*utaut_nsl[k])/(gplus+wplus));
        //      printf("ut_tot_ns=%f  ut_tot_sl=%f ut_tot_t=%f\n",utaut_nsl[k],utaut_sl[k],utaut_avg[k]);
        }

	printf("utau=%f\n",utauavg);

        uttotal = utaut_avg[up-2];
        printf("uttotal=%f error=%f\n",uttotal,fabs((utauavg-uttotal)/utauavg)*100.);
/*        for (k=0;k<(up-1);k++)
        {
                dsl[k] = (uttotal*uttotal - utaut_sl[k]*utaut_sl[k])*gplus;
                dnsl[k] = (utaut_nsl[k]*utaut_nsl[k] - uttotal*uttotal)*wplus;
                dtot[k] = 0.5*(dsl[k]+dnsl[k]);
        }
        for (k=0;k<(up-1);k++)
        {
                utaut_nsl[k] = sqrt(uttotal*uttotal + (dtot[k]/((double)wplus)));
                utaut_sl[k] = sqrt(uttotal*uttotal - (dtot[k]/((double)gplus)));
        //      printf("ut_tot_ns=%f  ut_tot_sl=%f ut_tot_t=%f\n",utaut_nsl[k],utaut_sl[k],utaut_avg[k]);
        }
*/
	
//********************************************************************************************		
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			if (S[j] == 0)
			{
				Trii1np[k] += (nu*Trii[j][k]/(wstravg[j]*wstravg[j]));
	                        Epii1np[k] += (nu*Epii[j][k]/(wstravg[j]*wstravg[j]));
	                        Pii1np[k] += (nu*Pii[j][k]/(wstravg[j]*wstravg[j]));
                        
                        
	                        P11np[k] += (nu*P1[j][k]/(wstravg[j]*wstravg[j]));
	                        P21np[k] += (nu*P2[j][k]/(wstravg[j]*wstravg[j]));
	                        P31np[k] += (nu*P3[j][k]/(wstravg[j]*wstravg[j]));
	                }
		}
	}

	for (k=0;k<zd;k++)
	{
		Trii1np[k] /= ((double)ct);
                Epii1np[k] /= ((double)ct);
                Pii1np[k] /= ((double)ct);
                
                P11np[k] /= ((double)ct);
                P21np[k] /= ((double)ct);
                P31np[k] /= ((double)ct);
	}
	printf("ct no-slip check %d\n",ct);

	sprintf(stri,TEXT);
	sprintf(fn,"Epii-no-slip-collected-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(Epii1np[k]+Epii1np[zd-1-k]));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"Pii-no-slip-collected-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pii1np+k))+(*(Pii1np+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"Trii-no-slip-collected-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Trii1np+k))+(*(Trii1np+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"P1-no-slip-collected-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P11np+k))+(*(P11np+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"P2-no-slip-collected-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P21np+k))+(*(P21np+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"P3-no-slip-collected-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P31np+k))+(*(P31np+zd-1-k))));
        fprintf(tcp,"\n");
        fclose(tcp);

//***************************************************************************************************
	sprintf(fn,"Epii-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(Epii1[k]+Epii1[zd-1-k])*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epii-combined-avg-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(Eps[k]+Eps[zd-1-k])*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epii-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),Epii1[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epii-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(Epii3[k]+Epii3[zd-1-k])*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epii-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),Epii3[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pii1+k))+(*(Pii1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Pii1+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pii3+k))+(*(Pii3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Pii3+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-combined-avg-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(Prd[k]+Prd[zd-1-k])*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trii-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Trii1+k))+(*(Trii1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tiir-no-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiir1+k))+(*(Tiir1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"Tiipai-no-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiipai1+k))+(*(Tiipai1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"Tiis-no-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiis1+k))+(*(Tiis1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"Trii-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Trii1+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trii-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Trii3+k))+(*(Trii3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tiir-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiir3+k))+(*(Tiir3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"Tiipai-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiipai3+k))+(*(Tiipai3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"Tiis-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiis3+k))+(*(Tiis3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"Trii-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Trii3+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(VdKdy1+k))+(*(VdKdy1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(VdKdy1+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(VdKdy3+k))+(*(VdKdy3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(VdKdy3+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(WdKdz1+k))+(*(WdKdz1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(WdKdz1+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(WdKdz3+k))+(*(WdKdz3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(WdKdz3+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P11+k))+(*(P11+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P11+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P13+k))+(*(P13+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P13+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P21+k))+(*(P21+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P21+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P23+k))+(*(P23+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P23+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P31+k))+(*(P31+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P31+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P33+k))+(*(P33+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P33+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P41+k))+(*(P41+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P41+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P43+k))+(*(P43+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P43+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P51+k))+(*(P51+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P51+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P53+k))+(*(P53+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P53+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P61+k))+(*(P61+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P61+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P63+k))+(*(P63+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P63+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"budget-slip-no-slip-collected-w-tau-total-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES =\"Z+\",\"Z+t\",\"Z/h\",\"Piin\",\"uvdUdyn\",\"uwdUdzn\""),
	fprintf(tcp,",\"Epiin\",\"Triin\"");
	fprintf(tcp,"\"Piis\",\"uvdUdys\",\"uwdUdzs\",\"Epiis\",\"Triis\"\n");
	for (k=1;k<(up-1);k++)
        {
                fprintf(tcp,"%f ",((*(z+k))+1)*Ret*utaut_sl[k]/utauavg);
                fprintf(tcp,"%f ",((*(z+k))+1)*Ret*utaut_nsl[k]/utauavg);
                fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%.9f ",0.5*((*(Pii1+k))+(*(Pii1+zd-1-k)))*nu/(utaut_nsl[k]*utaut_nsl[k]*utaut_nsl[k]*utaut_nsl[k]));
		fprintf(tcp,"%.9f ",0.5*((*(P11+k))+(*(P11+zd-1-k)))*nu/(utaut_nsl[k]*utaut_nsl[k]*utaut_nsl[k]*utaut_nsl[k]));
		fprintf(tcp,"%.9f ",0.5*((*(P21+k))+(*(P21+zd-1-k)))*nu/(utaut_nsl[k]*utaut_nsl[k]*utaut_nsl[k]*utaut_nsl[k]));
		fprintf(tcp,"%.9f ",0.5*(Epii1[k]+Epii1[zd-1-k])*nu/(utaut_nsl[k]*utaut_nsl[k]*utaut_nsl[k]*utaut_nsl[k]));
		fprintf(tcp,"%.9f ",0.5*((*(Trii1+k))+(*(Trii1+zd-1-k)))*nu/(utaut_nsl[k]*utaut_nsl[k]*utaut_nsl[k]*utaut_nsl[k]));
		fprintf(tcp,"%.9f ",0.5*((*(Pii3+k))+(*(Pii3+zd-1-k)))*nu/(utaut_sl[k]*utaut_sl[k]*utaut_sl[k]*utaut_sl[k]));
		fprintf(tcp,"%.9f ",0.5*((*(P13+k))+(*(P13+zd-1-k)))*nu/(utaut_sl[k]*utaut_sl[k]*utaut_sl[k]*utaut_sl[k]));
		fprintf(tcp,"%.9f ",0.5*((*(P23+k))+(*(P23+zd-1-k)))*nu/(utaut_sl[k]*utaut_sl[k]*utaut_sl[k]*utaut_sl[k]));
		fprintf(tcp,"%.9f ",0.5*(Epii3[k]+Epii3[zd-1-k])*nu/(utaut_sl[k]*utaut_sl[k]*utaut_sl[k]*utaut_sl[k]));
		fprintf(tcp,"%.9f\n",0.5*((*(Trii3+k))+(*(Trii3+zd-1-k)))*nu/(utaut_sl[k]*utaut_sl[k]*utaut_sl[k]*utaut_sl[k]));
	}
	fprintf(tcp,"\n");
        fclose(tcp);

	free(z);
}
