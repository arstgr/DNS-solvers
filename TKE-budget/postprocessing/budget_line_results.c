#include "definitions.h"

void budget_Line_results(int xd, int yd, int zd, int st, int en)
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

	DP Pii1[zd],Epii1[zd],Trii1[zd],VdKdy1[zd],WdKdz1[zd],U1[zd];
	DP Pii2[zd],Epii2[zd],Trii2[zd],VdKdy2[zd],WdKdz2[zd],U2[zd];
	DP Pii3[zd],Epii3[zd],Trii3[zd],VdKdy3[zd],WdKdz3[zd],U3[zd];
	
	DP P11[zd],P21[zd],P31[zd],P41[zd],P51[zd],P61[zd];
	DP P12[zd],P22[zd],P32[zd],P42[zd],P52[zd],P62[zd];
	DP P13[zd],P23[zd],P33[zd],P43[zd],P53[zd],P63[zd];

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
	printf("budget_Line_results: ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);

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

		Pii2[k] = 0.;
		Epii2[k] = 0.;
		Trii2[k] = 0.;
		VdKdy2[k] = 0.;
		WdKdz2[k] = 0.;
		U2[k] = 0.;
		P12[k] = 0.;
		P22[k] = 0.;
		P32[k] = 0.;
		P42[k] = 0.;
		P52[k] = 0.;
		P62[k] = 0.;

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
	}

	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			Pii[j][k] = 0.;
			Epii[j][a] = 0.;
			VdKdy[j][k] = 0.;
			WdKdz[j][k] = 0.;
			Trii[j][k] = 0.;

			P1[j][k] = 0.;
			P2[j][k] = 0.;
			P3[j][k] = 0.;
			P4[j][k] = 0.;
			P5[j][k] = 0.;
			P6[j][k] = 0.;
		}
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
			for (a=0;a<zd;a++)
				Trii[j][a] += Triit[j][a];
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(Temp,sizeof(DP),(yd*zd),sv);
		fread(&Triit[0][0],sizeof(DP),(yd*zd),sv);
		for (j=0;j<yd;j++)
			for (a=0;a<zd;a++)
				Trii[j][a] += Triit[j][a];
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
			for (a=0;a<zd;a++)
				Trii[j][a] += Triit[j][a];
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
		}
	}

	ct=0;
	for (j=0;j<yd;j+=(gplus+wplus))
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

			Trii2[k] += Trii[j+wst][k];
			Epii2[k] += Epii[j+wst][k];
			Pii2[k] += Pii[j+wst][k];
			VdKdy2[k] += VdKdy[j+wst][k];
			WdKdz2[k] += WdKdz[j+wst][k];
			P12[k] += P1[j+wst][k];
			P22[k] += P2[j+wst][k];
			P32[k] += P3[j+wst][k];
			P42[k] += P4[j+wst][k];
			P52[k] += P5[j+wst][k];
			P62[k] += P6[j+wst][k];

			Trii3[k] += Trii[j+wst+gplus/2][k];
			Epii3[k] += Epii[j+wst+gplus/2][k];
			Pii3[k] += Pii[j+wst+gplus/2][k];
			VdKdy3[k] += VdKdy[j+wst+gplus/2][k];
			WdKdz3[k] += WdKdz[j+wst+gplus/2][k];
			P13[k] += P1[j+wst+gplus/2][k];
			P23[k] += P2[j+wst+gplus/2][k];
			P33[k] += P3[j+wst+gplus/2][k];
			P43[k] += P4[j+wst+gplus/2][k];
			P53[k] += P5[j+wst+gplus/2][k];
			P63[k] += P6[j+wst+gplus/2][k];
		}
		ct++;
	}

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

		Trii2[k] /= ((double)ct);
		Epii2[k] /= ((double)ct);
		Pii2[k] /= ((double)ct);
		VdKdy2[k] /= ((double)ct);
		WdKdz2[k] /= ((double)ct);
		P12[k] /= ((double)ct);
		P22[k] /= ((double)ct);
		P32[k] /= ((double)ct);
		P42[k] /= ((double)ct);
		P52[k] /= ((double)ct);
		P62[k] /= ((double)ct);
	
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
	}

	sprintf(stri,TEXT);
	sprintf(fn,"Epii-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(Epii1[k]+Epii1[zd-1-k])*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epii-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),Epii1[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epii-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(Epii2[k]+Epii2[zd-1-k])*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epii-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),Epii2[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epii-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(Epii3[k]+Epii3[zd-1-k])*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epii-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),Epii3[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pii1+k))+(*(Pii1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Pii1+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pii2+k))+(*(Pii2+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Pii2+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pii3+k))+(*(Pii3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Pii3+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trii-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Trii1+k))+(*(Trii1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trii-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Trii1+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trii-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Trii2+k))+(*(Trii2+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trii-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Trii2+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trii-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Trii3+k))+(*(Trii3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trii-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Trii3+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(VdKdy1+k))+(*(VdKdy1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(VdKdy1+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(VdKdy2+k))+(*(VdKdy2+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(VdKdy2+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(VdKdy3+k))+(*(VdKdy3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(VdKdy3+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(WdKdz1+k))+(*(WdKdz1+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(WdKdz1+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(WdKdz2+k))+(*(WdKdz2+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(WdKdz2+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(WdKdz3+k))+(*(WdKdz3+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(WdKdz3+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P11+k))+(*(P11+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P11+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P12+k))+(*(P12+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P12+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P13+k))+(*(P13+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P13+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P21+k))+(*(P21+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P21+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P22+k))+(*(P22+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P22+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P23+k))+(*(P23+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P23+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P31+k))+(*(P31+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P31+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P32+k))+(*(P32+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P32+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P33+k))+(*(P33+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P33+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P41+k))+(*(P41+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P41+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P42+k))+(*(P42+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P42+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P43+k))+(*(P43+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P43+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P51+k))+(*(P51+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P51+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P52+k))+(*(P52+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P52+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P53+k))+(*(P53+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P53+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P61+k))+(*(P61+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P61+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P62+k))+(*(P62+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P62+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(P63+k))+(*(P63+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(P63+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	free(z);
}