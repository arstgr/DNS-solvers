#include "definitions.h"

void aniso_Line_results(int xd, int yd, int zd, int st, int en)
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
	DP utin[yd][zd],vtin[yd][zd],wtin[yd][zd],uw[yd][zd];
	DP utint[yd][zd],vtint[yd][zd],wtint[yd][zd],uwt[yd][zd];
	DP uv[yd][zd],uvt[yd][zd],vw[yd][zd],vwt[yd][zd];
	DP UW[yd][zd],UWt[yd][zd];
	DP Uavgt[yd][zd],Uavg[yd][zd];
	DP dUdz1[zd],dUdz2[zd],dUdz3[zd];
	DP du1[zd],du2[zd],du3[zd];
	DP U1[zd],U2[zd],U3[zd];
	DP KE,KEU;

	DP uu1[zd],uv1[zd],uw1[zd],vv1[zd],vw1[zd],ww1[zd],UW1[zd];
	DP uu2[zd],uv2[zd],uw2[zd],vv2[zd],vw2[zd],ww2[zd],UW2[zd];
	DP uu3[zd],uv3[zd],uw3[zd],vv3[zd],vw3[zd],ww3[zd],UW3[zd];

	DP A[3][3]={0.};
	DP IIa[zd],IIIa[zd];

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
	printf("aniso_Line_results: ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);

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
		uu1[k] = 0.;
		uv1[k] = 0.;
		uw1[k] = 0.;
		vv1[k] = 0.;
		vw1[k] = 0.;
		ww1[k] = 0.;
		UW1[k] = 0.;
		dUdz1[k] = 0.;
		U1[k] = 0.;

		uu2[k] = 0.;
		uv2[k] = 0.;
		uw2[k] = 0.;
		vv2[k] = 0.;
		vw2[k] = 0.;
		ww2[k] = 0.;
		UW2[k] = 0.;
		dUdz2[k] = 0.;
		U2[k] = 0.;

		uu3[k] = 0.;
		uv3[k] = 0.;
		uw3[k] = 0.;
		vv3[k] = 0.;
		vw3[k] = 0.;
		ww3[k] = 0.;
		UW3[k] = 0.;
		dUdz3[k] = 0.;
		U3[k] = 0.;
	}

	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			uv[j][k] = 0.;
			vw[j][a] = 0.;

			utin[j][k] = 0.;
			vtin[j][k] = 0.;
			wtin[j][k] = 0.;
		
			uw[j][k] = 0.;
			UW[j][k] = 0.;
			
			Uavg[j][k] = 0.;
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
		fread(&utint[0][0],sizeof(DP),(yd*zd),sv);
		fread(&vtint[0][0],sizeof(DP),(yd*zd),sv);
		fread(&wtint[0][0],sizeof(DP),(yd*zd),sv);
		fread(&uvt[0][0],sizeof(DP),(yd*zd),sv);
		fread(&uwt[0][0],sizeof(DP),(yd*zd),sv);
		fread(&vwt[0][0],sizeof(DP),(yd*zd),sv);
		fread(&UWt[0][0],sizeof(DP),(zd*yd),sv);
		fclose(sv);

		sprintf(fn,"uAVG-2d.%.3d.%.4d",0,k);
		sv=fopen(fn,"rb");
		fread(&Uavgt[0][0],sizeof(DP),yd*zd,sv);
		fclose(sv);

		for (j=0;j<yd;j++)
		{
			for (a=1;a<(zd-1);a++)
			{
				uv[j][a] += uvt[j][a];
				vw[j][a] += vwt[j][a];

				utin[j][a] += utint[j][a];
				vtin[j][a] += vtint[j][a];
				wtin[j][a] += wtint[j][a];
			
				uw[j][a] += uwt[j][a];
				UW[j][a] += UWt[j][a];

				Uavg[j][a] += Uavgt[j][a];
			}
		}
	}
	
	ct=0;
	ct=yd/(gplus+wplus);

	for (j=0;j<yd;j++)
	{	
		for (k=0;k<zd;k++)
		{
			uv[j][k] /= ((double)num);
			vw[j][k] /= ((double)num);
			UW[j][k] /= ((double)num);

			utin[j][k] /= ((double)num);
			vtin[j][k] /= ((double)num);
			wtin[j][k] /= ((double)num);
			uw[j][k] /= ((double)num);

			Uavg[j][k] /= ((double)num);
		}
	}

	for (j=0;j<yd;j+=(gplus+wplus))
	{
		for (k=0;k<zd;k++)
		{	
			uu1[k] += utin[j][k];
			uv1[k] += uv[j][k];
			uw1[k] += uw[j][k];
			vv1[k] += vtin[j][k];
			vw1[k] += vw[j][k];
			ww1[k] += wtin[j][k];
			UW1[k] += UW[j][k];
			U1[k] += Uavg[j][k];

			uu2[k] += utin[j+wst][k];
			uv2[k] += uv[j+wst][k];
			uw2[k] += uw[j+wst][k];
			vv2[k] += vtin[j+wst][k];
			vw2[k] += vw[j+wst][k];
			ww2[k] += wtin[j+wst][k];
			UW2[k] += UW[j+wst][k];
			U2[k] += Uavg[j+wst][k];

			uu3[k] += utin[j+wst+gplus/2][k];
			uv3[k] += uv[j+wst+gplus/2][k];
			uw3[k] += uw[j+wst+gplus/2][k];
			vv3[k] += vtin[j+wst+gplus/2][k];
			vw3[k] += vw[j+wst+gplus/2][k];
			ww3[k] += wtin[j+wst+gplus/2][k];
			UW3[k] += UW[j+wst+gplus/2][k];
			U3[k] += Uavg[j+wst+gplus/2][k];
		}
	}

	for (k=0;k<zd;k++)
	{	
		uu1[k] /= ((double)ct);
		uv1[k] /= ((double)ct);
		uw1[k] /= ((double)ct);
		vv1[k] /= ((double)ct);
		vw1[k] /= ((double)ct);
		ww1[k] /= ((double)ct);
		UW1[k] /= ((double)ct);
		U1[k] /= ((double)ct);

		uu2[k] /= ((double)ct);
		uv2[k] /= ((double)ct);
		uw2[k] /= ((double)ct);
		vv2[k] /= ((double)ct);
		vw2[k] /= ((double)ct);
		ww2[k] /= ((double)ct);
		UW2[k] /= ((double)ct);
		U2[k] /= ((double)ct);

		uu3[k] /= ((double)ct);
		uv3[k] /= ((double)ct);
		uw3[k] /= ((double)ct);
		vv3[k] /= ((double)ct);
		vw3[k] /= ((double)ct);
		ww3[k] /= ((double)ct);
		UW3[k] /= ((double)ct);
		U3[k] /= ((double)ct);
	}

	for (k=0;k<up;k++)
	{
                du1[k] = 0.5*(U1[k]+U1[zd-1-k]);
		du2[k] = 0.5*(U2[k]+U2[zd-1-k]);
		du3[k] = 0.5*(U3[k]+U3[zd-1-k]);
	}

        dUdz1[0] = (-2.*du1[1]+3.*du1[2]-du1[3]);
	dUdz2[0] = (-2.*du2[1]+3.*du2[2]-du2[3]);
	dUdz3[0] = (-2.*du3[1]+3.*du3[2]-du3[3]);
        
        dUdz1[1] = (-3.*du1[1]+4.*du1[2]-du1[3])/2.;
	dUdz2[1] = (-3.*du2[1]+4.*du2[2]-du2[3])/2.;
	dUdz3[1] = (-3.*du3[1]+4.*du3[2]-du3[3])/2.;

        for (k=2;k<(up-2);k++)
	{
                dUdz1[k] = (-3.*du1[k]+4.*du1[k+1]-du1[k+2])/2.;
		dUdz2[k] = (-3.*du2[k]+4.*du2[k+1]-du2[k+2])/2.;
		dUdz3[k] = (-3.*du3[k]+4.*du3[k+1]-du3[k+2])/2.;
	}
        dUdz1[up-2]=du1[up-2]-du1[up-3];
        dUdz1[up-1]=du1[up-1]-du1[up-2];

	dUdz2[up-2]=du2[up-2]-du2[up-3];
        dUdz2[up-1]=du2[up-1]-du2[up-2];

	dUdz3[up-2]=du3[up-2]-du3[up-3];
        dUdz3[up-1]=du3[up-1]-du3[up-2];
			

	for (k=0;k<zd;k++)
	{
		IIa[k] = 0.;
		IIIa[k] = 0.;
	}

	for (k=1;k<up;k++)
	{
		KE = (uu1[k]+vv1[k]+ww1[k]);
		KEU = uu1[zd-1-k]+vv1[zd-1-k]+ww1[zd-1-k];
		A[0][0] = 0.5*((uu1[k]/KE)+(uu1[zd-1-k]/KEU))-1./3.;
		A[0][1] = 0.5*((uv1[k]/KE)-(uv1[zd-1-k]/KEU));
		A[1][0] = A[0][1];
		A[0][2] = 0.5*((uw1[k]/KE)-(uw1[zd-1-k]/KEU));
		A[2][0] = A[0][2];
		A[1][1] = 0.5*((vv1[k]/KE)+(vv1[zd-1-k]/KEU))-1./3.;
		A[1][2] = 0.5*((vw1[k]/KE)-(vw1[zd-1-k]/KEU));
		A[2][1] = A[1][2];
		A[2][2] = 0.5*((ww1[k]/KE)+(ww1[zd-1-k]/KEU))-1./3.;

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

	sprintf(fn,"Aniso-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"II\",\"III\"\n");
	for (k=up-1;k>0;k--)
		fprintf(tcp,"%f %.9f %.9f\n",((*(z+k))+1)*Ret,IIa[k],IIIa[k]);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-no-slip-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"\n");
//	fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
	for (k=0;k<up;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",0.5*(U1[k]+U1[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(uu1[k]+uu1[zd-1-k]))/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(vv1[k]+vv1[zd-1-k]))/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(ww1[k]+ww1[zd-1-k]))/utauavg);
		fprintf(tcp,"%f ",0.5*(uv1[k]-uv1[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(uw1[k]-uw1[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(vw1[k]-vw1[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",0.5*(UW1[k]-UW1[zd-1-k])/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-no-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"\n");
//	fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
	for (k=0;k<zd;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",(U1[k])/utauavg);
		fprintf(tcp,"%f ",sqrt(uu1[k])/utauavg);
		fprintf(tcp,"%f ",sqrt(vv1[k])/utauavg);
		fprintf(tcp,"%f ",sqrt(ww1[k])/utauavg);
		fprintf(tcp,"%f ",(uv1[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(uw1[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(vw1[k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",(UW1[k])/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	for (k=0;k<zd;k++)
	{
		IIa[k] = 0.;
		IIIa[k] = 0.;
	}

	for (k=1;k<up;k++)
	{
		KE = (uu2[k]+vv2[k]+ww2[k]);
		KEU = uu2[zd-1-k]+vv2[zd-1-k]+ww2[zd-1-k];
		A[0][0] = 0.5*((uu2[k]/KE)+(uu2[zd-1-k]/KEU))-1./3.;
		A[0][1] = 0.5*((uv2[k]/KE)-(uv2[zd-1-k]/KEU));
		A[1][0] = A[0][1];
		A[0][2] = 0.5*((uw2[k]/KE)-(uw2[zd-1-k]/KEU));
		A[2][0] = A[0][2];
		A[1][1] = 0.5*((vv2[k]/KE)+(vv2[zd-1-k]/KEU))-1./3.;
		A[1][2] = 0.5*((vw2[k]/KE)-(vw2[zd-1-k]/KEU));
		A[2][1] = A[1][2];
		A[2][2] = 0.5*((ww2[k]/KE)+(ww2[zd-1-k]/KEU))-1./3.;

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

	sprintf(fn,"Aniso-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"II\",\"III\"\n");
	for (k=up-1;k>0;k--)
		fprintf(tcp,"%f %.9f %.9f\n",((*(z+k))+1)*Ret,IIa[k],IIIa[k]);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-edge-slip-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"\n");
//	fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
	for (k=0;k<up;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",0.5*(U2[k]+U2[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(uu2[k]+uu2[zd-1-k]))/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(vv2[k]+vv2[zd-1-k]))/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(ww2[k]+ww2[zd-1-k]))/utauavg);
		fprintf(tcp,"%f ",0.5*(uv2[k]-uv2[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(uw2[k]-uw2[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(vw2[k]-vw2[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",0.5*(UW2[k]-UW2[zd-1-k])/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-edge-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"\n");
//	fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
	for (k=0;k<zd;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",(U2[k])/utauavg);
		fprintf(tcp,"%f ",sqrt(uu2[k])/utauavg);
		fprintf(tcp,"%f ",sqrt(vv2[k])/utauavg);
		fprintf(tcp,"%f ",sqrt(ww2[k])/utauavg);
		fprintf(tcp,"%f ",(uv2[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(uw2[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(vw2[k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",(UW2[k])/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	for (k=0;k<zd;k++)
	{
		IIa[k] = 0.;
		IIIa[k] = 0.;
	}

	for (k=1;k<up;k++)
	{
		KE = (uu3[k]+vv3[k]+ww3[k]);
		KEU = uu3[zd-1-k]+vv3[zd-1-k]+ww3[zd-1-k];
		A[0][0] = 0.5*((uu3[k]/KE)+(uu3[zd-1-k]/KEU))-1./3.;
		A[0][1] = 0.5*((uv3[k]/KE)-(uv3[zd-1-k]/KEU));
		A[1][0] = A[0][1];
		A[0][2] = 0.5*((uw3[k]/KE)-(uw3[zd-1-k]/KEU));
		A[2][0] = A[0][2];
		A[1][1] = 0.5*((vv3[k]/KE)+(vv3[zd-1-k]/KEU))-1./3.;
		A[1][2] = 0.5*((vw3[k]/KE)-(vw3[zd-1-k]/KEU));
		A[2][1] = A[1][2];
		A[2][2] = 0.5*((ww3[k]/KE)+(ww3[zd-1-k]/KEU))-1./3.;

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

	sprintf(fn,"Aniso-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"II\",\"III\"\n");
	for (k=up-1;k>0;k--)
		fprintf(tcp,"%f %.9f %.9f\n",((*(z+k))+1)*Ret,IIa[k],IIIa[k]);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-slip-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"\n");
//	fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
	for (k=0;k<up;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",0.5*(U3[k]+U3[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(uu3[k]+uu3[zd-1-k]))/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(vv3[k]+vv3[zd-1-k]))/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(ww3[k]+ww3[zd-1-k]))/utauavg);
		fprintf(tcp,"%f ",0.5*(uv3[k]-uv3[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(uw3[k]-uw3[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(vw3[k]-vw3[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",0.5*(UW3[k]-UW3[zd-1-k])/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-slip-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"\n");
//	fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
	for (k=0;k<zd;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",(U3[k])/utauavg);
		fprintf(tcp,"%f ",sqrt(uu3[k])/utauavg);
		fprintf(tcp,"%f ",sqrt(vv3[k])/utauavg);
		fprintf(tcp,"%f ",sqrt(ww3[k])/utauavg);
		fprintf(tcp,"%f ",(uv3[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(uw3[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(vw3[k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",(UW3[k])/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	free(z);
}