#include "definitions.h"

void aniso_slip_noslip_results(int xd, int yd, int zd, int st, int en)
{
	int i,j,k,a;
	FILE *sv;
	FILE *tcp,*wstr;
	char fn[150];
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
	DP UA[yd][zd],UAt[yd][zd];
	DP U1[zd],U2[zd];
	DP dAVG1[zd],dAVG2[zd],dAVG[yd][zd],dAVGt[yd][zd];
	DP dfluc1[zd],dfluc2[zd],dfluc[yd][zd],dfluct[yd][zd];
	DP dflucsq1[zd],dflucsq2[zd],dflucsq[yd][zd],dflucsq3[zd];
	DP KE,KEU;
	DP wstrs[yd],u1,u2,u3,utns=0.,u4,u5,up1,tmp;
	DP dusl[zd],dunsl[zd],dudxsl[zd],dudxnsl[zd];
	DP utinsq[yd][zd],vtinsq[yd][zd],wtinsq[yd][zd];
	double ubslip=0., ubnslip=0., ubtemp=0.,uavg[zd],rmstesta,rmstestb,rmstestc;
	double absuvsl[zd],absuvnsl[zd],utsl=0.,uttot=0.;
	double utaut_nsl[zd],utaut_sl[zd],utaut_avg[zd]; 
	DP uttotal,dsl[zd],dnsl[zd],dtot[zd];
	DP dutdz[zd],function[zd],integ[zd];

	DP uu1[zd],uv1[zd],uw1[zd],vv1[zd],vw1[zd],ww1[zd],UW1[zd],uusq1[zd],vvsq1[zd],wwsq1[zd];
	DP uu2[zd],uv2[zd],uw2[zd],vv2[zd],vw2[zd],ww2[zd],UW2[zd],uusq2[zd],vvsq2[zd],wwsq2[zd];
	DP uusq3[zd],vvsq3[zd],wwsq3[zd],uu3[zd],vv3[zd],ww3[zd],UW3[zd];	

	DP u_a,utin_a,vtin_a,wtin_a,uw_a,UW_a,uv_a,vw_a;
	DP U1np[zd],uu1np[zd],vv1np[zd],ww1np[zd],uv1np[zd],uw1np[zd],vw1np[zd],UW1np[zd],dAVG1np[zd],dfluc1np[zd];
	int SS[yd];
	DP avgdensity=0.;	

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
	printf("aniso_slip_noslip_results: ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);

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
		U1[k] = 0.;

		uu1np[k] = 0.;
                uv1np[k] = 0.;
                uw1np[k] = 0.;
                vv1np[k] = 0.;
                vw1np[k] = 0.;
                ww1np[k] = 0.;
                UW1np[k] = 0.;
                U1np[k] = 0.;

		uu2[k] = 0.;
		uv2[k] = 0.;
		uw2[k] = 0.;
		vv2[k] = 0.;
		vw2[k] = 0.;
		ww2[k] = 0.;
		UW2[k] = 0.;
		U2[k] = 0.;

		uu3[k] = 0.;
                vv3[k] = 0.;
                ww3[k] = 0.;

		UW3[k] = 0.;

		uusq1[k] = 0.;
		vvsq1[k] = 0.;
		wwsq1[k] = 0.;

		uusq2[k] = 0.;
		vvsq2[k] = 0.;
		wwsq2[k] = 0.;

		uusq3[k] = 0.;
		vvsq3[k] = 0.;
		wwsq3[k] = 0.;

		dusl[k] = 0.;
		dunsl[k] = 0.;
		dudxsl[k] = 0.;
		dudxnsl[k] = 0.;

		dAVG1[k] = 0.;
		dAVG1np[k] = 0.;
		dAVG2[k] = 0.;

		dfluc1[k] = 0.;
		dfluc1np[k] = 0.;
		dfluc2[k] = 0.;

		dflucsq1[k] = 0.;
		dflucsq2[k] = 0.;
		dflucsq3[k] = 0.;

		uavg[k] = 0.;

		absuvsl[k] = 0.;
		absuvnsl[k] = 0.;
	}

	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			uv[j][k] = 0.;
			vw[j][k] = 0.;

			utin[j][k] = 0.;
			vtin[j][k] = 0.;
			wtin[j][k] = 0.;

			utinsq[j][k] = 0.;
			vtinsq[j][k] = 0.;
			wtinsq[j][k] = 0.;
		
			uw[j][k] = 0.;
			UW[j][k] = 0.;

			UA[j][k] = 0.;
			UAt[j][k] = 0.;

			dAVG[j][k] = 0.;
			dfluc[j][k] = 0.;
			dflucsq[j][k] = 0.;

			dAVGt[j][k] = 0.;
                        dfluct[j][k] = 0.;
		}
	}
	
	for (j=0;j<yd;j++)
	{
		SS[j] = 5;
		wstrs[j] = 0.;
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
		fread(&dfluct[0][0],sizeof(DP),(zd*yd),sv);
		fread(&dAVGt[0][0],sizeof(DP),yd*zd,sv);
		fclose(sv);

		sprintf(fn,"uAVG-2d.%.3d.%.4d",0,k);
		sv=fopen(fn,"rb");
		fread(&UAt[0][0],sizeof(DP),yd*zd,sv);
		fclose(sv);
/*
		sprintf(fn,"dAVG-2d.%.3d.%.4d",0,k);
                sv=fopen(fn,"rb");
                fread(&dAVGt[0][0],sizeof(DP),yd*zd,sv);
                fclose(sv);

		sprintf(fn,"dens-fluc-2d.%.3d.%.4d",0,k);
                sv=fopen(fn,"rb");
                fread(&dfluct[0][0],sizeof(DP),yd*zd,sv);
                fclose(sv);
*/
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

				UA[j][a] += UAt[j][a];

				dfluc[j][a] += dfluct[j][a];
				dAVG[j][a] += dAVGt[j][a];
			}
		}
/*		ubtemp=0.;
		for (a=0;a<zd;a++)
			for (j=0;j<yd;j++)
				ubtemp += UAt[j][a];
		printf("time=%d ubtemp=%f ub=%f\n",k,ubtemp/((double)(yd*(zd-2.))),ub);
*/
	}

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

			UA[j][k] /= ((double)num);

			dfluc[j][k] /= ((double)num);
                        dAVG[j][k] /= ((double)num);
		}
	}

	for (j=0;j<yd;j++)
		for (k=1;k<(zd-1);k++)
			avgdensity += dAVG[j][k];

	avgdensity /= ((double)(yd*(zd-2.)));
	printf("average density =%f \n",avgdensity);

	for (j=0;j<yd;j++)
	{	
		for (k=0;k<zd;k++)
		{
			utinsq[j][k] = sqrt(utin[j][k]);
			vtinsq[j][k] = sqrt(vtin[j][k]);
			wtinsq[j][k] = sqrt(wtin[j][k]);
			dflucsq[j][k] = sqrt(dfluc[j][k]);
		}
	}

	for (k=0;k<zd;k++)
	{
		for (j=0;j<yd;j++)
		{
			uavg[k] += UA[j][k];
		}
		uavg[k] /= ((double)yd);
	}

	for (j=0;j<yd;j++)
	{	
		for (k=0;k<zd;k++)
		{
			uusq3[k] += utinsq[j][k];
			vvsq3[k] += vtinsq[j][k];
			wwsq3[k] += wtinsq[j][k];
			dflucsq3[k] += dflucsq[j][k];

			uu3[k] += utin[j][k];
			vv3[k] += vtin[j][k];
			ww3[k] += wtin[j][k];

			UW3[k] += UW[j][k];
		}
	}

	for (k=0;k<zd;k++)
	{
		uusq3[k] /= ((double)yd);
		vvsq3[k] /= ((double)yd);
		wwsq3[k] /= ((double)yd);
		dflucsq3[k] /= ((double)yd);

		uu3[k] /= ((double)yd);
		vv3[k] /= ((double)yd);
		ww3[k] /= ((double)yd);

		UW3[k] /= ((double)yd);
	}


	for (j=0;j<yd;j++)
	{
		u1=0.5*(UA[j][1]+UA[j][zd-1-1]);
		u2=0.5*(UA[j][2]+UA[j][zd-1-2]);
		u3=0.5*(UA[j][3]+UA[j][zd-1-3]);
		u4=0.5*(UA[j][4]+UA[j][zd-1-4]);
		u5=0.5*(UA[j][5]+UA[j][zd-1-5]);

		wstrs[j] = nu*(-2.*u1+3.*u2-u3);
//		wstrs[j] = nu*(3.*u1-u2/3.);
//		wstrs[j] = nu*2.*u1;
//		tmp = (2./3.)*(u4-u2)-(1./12.)*(u5-u1);
//		up1 = -(7./3.)*u1+6.*u2-3.*u3-(2./3.)*u4+3.*tmp;
//		wstrs[j] = nu*((9./2.)*u1+(1./6.)*u2-(3./2.)*up1);
	}
	for (j=wst;j<yd;j+=(gplus+wplus))
		for (a=0;a<gplus;a++)
			wstrs[j+a] = 0.;

	for (j=0,uttot=0.;j<yd;j++)
		uttot+= wstrs[j];

	uttot /= ((double)yd);	

	sprintf(fn,"w-str-spnws-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Y+\",\"W\"\n");
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",(j+0.5)*utauavg/nu,wstrs[j]/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	for (j=wst;j<yd;j+=(gplus+wplus))
		for (a=0;a<gplus;a++)
			SS[j+a] = 0;

	ct=0;
	for (j=wst;j<yd;j+=(gplus+wplus))
	{
		for (a=0;a<gplus;a++)
		{
			for (k=0;k<zd;k++)
			{
				uu2[k] += utin[j+a][k];
				uv2[k] += uv[j+a][k];
				uw2[k] += uw[j+a][k];
				vv2[k] += vtin[j+a][k];
				vw2[k] += vw[j+a][k];
				ww2[k] += wtin[j+a][k];
				UW2[k] += UW[j+a][k];

				U2[k] += UA[j+a][k];

				dAVG2[k] += dAVG[j+a][k];
				dfluc2[k] += dfluc[j+a][k];

				uusq2[k] += utinsq[j+a][k];
				vvsq2[k] += vtinsq[j+a][k];
				wwsq2[k] += wtinsq[j+a][k];
				dflucsq2[k] += dflucsq[j+a][k];

				absuvsl[k] += fabs(uv[j+a][k]);
			}
			utsl+=wstrs[j+a];
			ct++;
		}
	}
	for (k=0;k<zd;k++)
		ubtemp += U2[k];
	for (k=0;k<zd;k++)
	{	
		uu2[k] /= ((double)ct);
		uv2[k] /= ((double)ct);
		uw2[k] /= ((double)ct);
		vv2[k] /= ((double)ct);
		vw2[k] /= ((double)ct);
		ww2[k] /= ((double)ct);
		UW2[k] /= ((double)ct);
		U2[k] /= ((double)ct);

		dAVG2[k] /= ((double)ct);
                dfluc2[k] /= ((double)ct);

		uusq2[k] /= ((double)ct);
		vvsq2[k] /= ((double)ct);
		wwsq2[k] /= ((double)ct);
		dflucsq2[k] /= ((double)ct);
		
		absuvsl[k] /= ((double)ct);
	}
	utsl /= ((double)ct);
	printf("ct slip: %d\n",ct);

	ct=0;
	for (j=(gplus+wst);j<(yd-wst);j+=(gplus+wplus))
	{
		for (a=0;a<wplus;a++)
		{
			for (k=0;k<zd;k++)
			{
				uu1[k] += utin[j+a][k];
				uv1[k] += uv[j+a][k];
				uw1[k] += uw[j+a][k];
				vv1[k] += vtin[j+a][k];
				vw1[k] += vw[j+a][k];
				ww1[k] += wtin[j+a][k];
				UW1[k] += UW[j+a][k];
	
				U1[k] += UA[j+a][k];

				dAVG1[k] += dAVG[j+a][k];
				dfluc1[k] += dfluc[j+a][k];

				uusq1[k] += utinsq[j+a][k];
				vvsq1[k] += vtinsq[j+a][k];
				wwsq1[k] += wtinsq[j+a][k];
				dflucsq1[k] += dflucsq[j+a][k];

				absuvnsl[k] += fabs(uv[j+a][k]);
			}
			utns += wstrs[j+a];
			ct++;
		}
	}
	for (a=0;a<wst;a++)
	{
		for (k=0;k<zd;k++)
		{
			uu1[k] += utin[a][k];
			uv1[k] += uv[a][k];
			uw1[k] += uw[a][k];
			vv1[k] += vtin[a][k];
			vw1[k] += vw[a][k];
			ww1[k] += wtin[a][k];
			UW1[k] += UW[a][k];
	
			U1[k] += UA[a][k];
			absuvnsl[k] += fabs(uv[a][k]);

			dAVG1[k] += dAVG[a][k];
                        dfluc1[k] += dfluc[a][k];

                        uusq1[k] += utinsq[a][k];
                        vvsq1[k] += vtinsq[a][k];
                        wwsq1[k] += wtinsq[a][k];
                        dflucsq1[k] += dflucsq[a][k];

			uu1[k] += utin[yd-wst+a][k];
			uv1[k] += uv[yd-wst+a][k];
			uw1[k] += uw[yd-wst+a][k];
			vv1[k] += vtin[yd-wst+a][k];
			vw1[k] += vw[yd-wst+a][k];
			ww1[k] += wtin[yd-wst+a][k];
			UW1[k] += UW[yd-wst+a][k];
	
			U1[k] += UA[yd-wst+a][k];

			dAVG1[k] += dAVG[yd-wst+a][k];
                        dfluc1[k] += dfluc[yd-wst+a][k];
			
			uusq1[k] += utinsq[yd-wst+a][k];
			vvsq1[k] += vtinsq[yd-wst+a][k];
			wwsq1[k] += wtinsq[yd-wst+a][k];
			dflucsq1[k] += dflucsq[yd-wst+a][k];
			
			absuvnsl[k] += fabs(uv[yd-wst+a][k]);
		}
		utns += wstrs[a];
		utns += wstrs[yd-wst+a];
		ct +=2;
	}

	for (k=0;k<zd;k++)
                ubtemp += U1[k];

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

		dAVG1[k] /= ((double)ct);
                dfluc1[k] /= ((double)ct);

		uusq1[k] /= ((double)ct);
		vvsq1[k] /= ((double)ct);
		wwsq1[k] /= ((double)ct);
		dflucsq1[k] /= ((double)ct);

		absuvnsl[k] /= ((double)ct);
	}
	utns /= ((double)ct);
	printf("utau-no-slip=%f, tau_ns_w=%.10f tau_sl_w=%.10f, ut_sl=%f, ut/ut-ns=%f utns/ut=%f   ubtemp=%f\n",sqrt(utns),utns,utsl,sqrt(utsl),utauavg/sqrt(utns),sqrt(utns)/utauavg,ubtemp/((double)(yd*(zd-2.))));
	printf("ct no-slip: %d\n",ct);
	printf("ut_w_calc=%.12f uttot=%.12f\n",sqrt((gplus*utsl+wplus*utns)/(gplus+wplus)),sqrt(uttot));
	ubtemp =0.;
	rmstesta=0.;
	rmstestb=0.;
	rmstestc=0.;
	for (k=0;k<zd;k++)
	{
		ubslip += U2[k];
		ubnslip += U1[k];
		ubtemp += uavg[k];

		rmstesta += uusq1[k];
		rmstestb += uusq2[k];
		rmstestc += uusq3[k];
	}
	ubslip /= ((double)(zd-2.));
	ubnslip /= ((double)(zd-2.));
	ubtemp /= ((double)(zd-2.));

	rmstesta /= ((double)(zd-2.));
	rmstestb /= ((double)(zd-2.));
	rmstestc /= ((double)(zd-2.));

	printf("ub-slip= %f       ub-noslip= %f     should be %f  ub=%f     ubtemp=%f\n",ubslip,ubnslip,(gplus*ubslip+wplus*ubnslip)/(gplus+wplus),ub,ubtemp);
	printf("rms-sq-first no-slip:%f slip:%f total %f should be %f\n",rmstesta,rmstestb,rmstestc,(gplus*rmstestb+wplus*rmstesta)/(gplus+wplus));

	rmstesta=0.;
        rmstestb=0.;
        rmstestc=0.;
	for (k=0;k<zd;k++)
        {
		rmstesta += uu1[k];
                rmstestb += uu2[k];
                rmstestc += uu3[k];
	}
	rmstesta /= ((double)(zd-2.));
        rmstestb /= ((double)(zd-2.));
        rmstestc /= ((double)(zd-2.));
	printf("rms-regular no-slip:%f slip:%f total %f should be %f\n",sqrt(rmstesta),sqrt(rmstestb),sqrt(rmstestc),sqrt((gplus*rmstestb+wplus*rmstesta)/(gplus+wplus)));

	for (k=0;k<up+1;k++)
	{
                dusl[k] = 0.5*(U2[k]+U2[zd-1-k]);
		dunsl[k] = 0.5*(U1[k]+U1[zd-1-k]);
	}

	for (k=0;k<up;k++)
	{
		dusl[zd-1-k] = dusl[k];
                dunsl[zd-1-k] = dunsl[k];
	}
/*
        dudxsl[0] = (-2.*dusl[1]+3.*dusl[2]-dusl[3]);
	dudxnsl[0] = (-2.*dunsl[1]+3.*dunsl[2]-dunsl[3]);
       
        dudxsl[1] = (-3.*dusl[1]+4.*dusl[2]-dusl[3])/2.;
	dudxnsl[1] = (-3.*dunsl[1]+4.*dunsl[2]-dunsl[3])/2.;

        for (k=2;k<(up-2);k++)
        {
		dudxsl[k] = (-3.*dusl[k]+4.*dusl[k+1]-dusl[k+2])/2.;
		dudxnsl[k] = (-3.*dunsl[k]+4.*dunsl[k+1]-dunsl[k+2])/2.;
	}

        dudxsl[up-2] = (dusl[up-2]-dusl[up-3]);
        dudxsl[up-1] = (dusl[up-1]-dusl[up-2]);

	dudxnsl[up-2] = (dunsl[up-2]-dunsl[up-3]);
        dudxnsl[up-1] = (dunsl[up-1]-dunsl[up-2]);
*/
// derivative based on 4th order compact FD scheme
	for (k=3;k<up;k++)
	{
		dudxsl[k] = (2./3.)*(dusl[k+1]-dusl[k-1])-(1./12.)*(dusl[k+2]-dusl[k-2]);
                dudxnsl[k] = (2./3.)*(dunsl[k+1]-dunsl[k-1])-(1./12.)*(dunsl[k+2]-dunsl[k-2]);
	}

	dudxsl[0] = 0.;
	dudxnsl[0] = 0.;

	dudxnsl[0] = 2.*dunsl[1];
	
//	dudxnsl[1] = 0.5*(-3.*dunsl[1]+4.*dunsl[2]-dunsl[3]);
//	dudxsl[1] = 0.5*(-3.*dusl[1]+4.*dusl[2]-dusl[3]);
	dudxsl[1] = 0.5*(dusl[2]-dusl[1]);
	dudxnsl[1] = 0.5*(dunsl[2]+dunsl[1]);
	
	dudxsl[2] = 0.5*(dusl[3]-dusl[1]);
	dudxnsl[2] = 0.5*(dunsl[3]-dunsl[1]);
	
//	dudxsl[1] = (2./3.)*(dusl[2]-dusl[1]) - (1./12.)*(dusl[3]-dusl[2]);
//	dudxnsl[1] = (2./3.)*(dunsl[2]+dunsl[1]) - (1./12.)*(dunsl[3]+dunsl[2]);

//	dudxsl[2] = 0.75*(dusl[3]-dusl[1]) - 0.25*dudxsl[1] - 0.25*dudxsl[3];
//	dudxnsl[2] = 0.75*(dunsl[3]-dunsl[1]) - 0.25*dudxnsl[1] - 0.25*dudxnsl[3];

	sprintf(fn,"dudx-slip-no-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
//        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"z<sup>+</sup>\", \"z/h\", \"<<greek>t</greek><sub>v</sub>><sup>+</sup><sub>ns</sub>\", \"<<greek>t</greek><sub>v</sub>><sup>+</sup><sub>sl</sub>\"");
	fprintf(tcp,",\"<<greek>t</greek><sub>v</sub>><sup>+<greek>*</greek></sup><sub>ns</sub>\",\"<<greek>t</greek><sub>v</sub>><sup>+<greek>*</greek></sup><sub>sl</sub>\",\"z<sup>+<greek>*</greek></sup>\"\n");
//        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
        {        
		fprintf(tcp,"%f %f %.9f %.9f",((*(z+k))+1)*Ret,((*(z+k))+1),nu*dudxnsl[k]/(utauavg*utauavg),nu*dudxsl[k]/(utauavg*utauavg));
		fprintf(tcp," %.9f %.9f %f\n",nu*dudxnsl[k]/(utauavg*utauavg*((gplus+wplus)/wplus)),nu*dudxsl[k]/(utauavg*utauavg*((gplus+wplus)/wplus)),((*(z+k))+1)*Ret*sqrt((double)(gplus+wplus)/wplus));
	}
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"stress-slip-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES = \"z<sup>+</sup>\", \"z/h\", \"t<sub>v,sl</sub><sup>+</sup>\", \"t<sub>R,sl</sub><sup>+</sup>\", \"t<sub>c,sl</sub><sup>+</sup>\"");
	fprintf(tcp," ,\"t<sub>v,nsl</sub><sup>+</sup>\", \"t<sub>R,nsl</sub><sup>+</sup>\", \"t<sub>c,nsl</sub><sup>+</sup>\"\n");
	for (k=0;k<up;k++)
        {
		fprintf(tcp,"%f %f %.9f %.9f %.9f",((*(z+k))+1)*Ret,((*(z+k))+1),nu*dudxsl[k]/(utauavg*utauavg),0.5*(uw2[k]-uw2[zd-1-k])/(utauavg*utauavg),0.5*(UW2[k]-UW2[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp," %.9f %.9f %.9f\n",nu*dudxnsl[k]/(utauavg*utauavg),0.5*(uw1[k]-uw1[zd-1-k])/(utauavg*utauavg),0.5*(UW1[k]-UW1[zd-1-k])/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
        fclose(tcp);

	for (k=0;k<(up-1);k++)
	{
		utaut_nsl[k] = sqrt((nu*dudxnsl[k]-0.5*(uw1[k]-uw1[zd-1-k]) )/(1.-((*(z+k))+1)));
		utaut_sl[k] = sqrt((nu*dudxsl[k] - 0.5*(uw2[k]-uw2[zd-1-k]) )/(1.-((*(z+k))+1)));
		utaut_avg[k] = sqrt((gplus*utaut_sl[k]*utaut_sl[k]+wplus*utaut_nsl[k]*utaut_nsl[k])/(gplus+wplus));
		printf("ut_tot_ns=%f  ut_tot_sl=%f ut_tot_t=%f\n",utaut_nsl[k],utaut_sl[k],utaut_avg[k]);
	}

	for (k=0;k<up;k++)
	{
		printf("ut_mod_nsl=%f ut_ns_org=%f ", sqrt(fabs(nu*dudxnsl[0]-(nu*dudxnsl[k]-0.5*(uw1[k]-uw1[zd-1-k])))/((*(z+k))+1)),utaut_nsl[k]);
		printf("ut_mod_sl =%f ut_sl_org=%f\n",sqrt(fabs(nu*dudxsl[0]-(nu*dudxsl[k]-0.5*(uw2[k]-uw2[zd-1-k])))/((*(z+k))+1)),utaut_sl[k]);
	}

	printf("utau=%f\n",utauavg); 
	
	uttotal = utaut_avg[up-2];

	dutdz[0] = 2.*utaut_nsl[1];
	dutdz[1] = 0.5*(utaut_nsl[2]+utaut_nsl[1]);
	dutdz[2] = 0.5*(utaut_nsl[3]-utaut_nsl[1]);
	for(k=3;k<(up-3);k++)
		dutdz[k] = (2./3.)*(utaut_nsl[k+1]-utaut_nsl[k-1])-(1./12.)*(utaut_nsl[k+2]-utaut_nsl[k-2]);
	dutdz[up-3] = 0.5*(utaut_nsl[up-2]-utaut_nsl[up-4]);
	dutdz[up-2] = utaut_nsl[up-2]-utaut_nsl[up-3];

	for (k=0;k<(up-1);k++)
		function[k] = log(k-0.5)*dutdz[k];

	function[0] = 0.;
	for (k=1;k<(up-1);k++)
		integ[k] = integral(function,k);

	integ[0] = 0.;
	
/*	uttotal = utauavg;
	printf("uttotal=%f error=%f\n",uttotal,fabs((utauavg-uttotal)/utauavg)*100.);
	for (k=0;k<(up-1);k++)
	{
		dsl[k] = (uttotal*uttotal - utaut_sl[k]*utaut_sl[k])*gplus + gplus*UW3[k];
		dnsl[k] = (utaut_nsl[k]*utaut_nsl[k] - uttotal*uttotal)*wplus - wplus*UW3[k];
		dtot[k] = 0.5*(dsl[k]+dnsl[k]);
	}
	for (k=0;k<(up-1);k++)
	{
		utaut_nsl[k] = sqrt(uttotal*uttotal + (dtot[k]/((double)wplus)));
		utaut_sl[k] = sqrt(uttotal*uttotal - (dtot[k]/((double)gplus)));
	//	printf("ut_tot_ns=%f  ut_tot_sl=%f ut_tot_t=%f\n",utaut_nsl[k],utaut_sl[k],utaut_avg[k]);
	}
*/
/*
	for (k=0;k<(up-1);k++)
	{
		utaut_nsl[k] = sqrt(utaut_nsl[k]*utaut_nsl[k] - ((double)(gplus+wplus))/((double)wplus) * UW3[k]);
		utaut_sl[k] = sqrt(utaut_sl[k]*utaut_sl[k] - ((double)(gplus+wplus))/((double)gplus) * UW3[k]);
	}
*/		
	sprintf(fn,"stress-slip-no-slip-mod-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"z<sup>+</sup>\", \"z/h\",\"t<sub>t,sl</sub><sup>+</sup>\", \"t<sub>t,nsl</sub><sup>+</sup>\"\n");
        for (k=0;k<(up-1);k++)
        {
                fprintf(tcp,"%f %f %.9f %.9f\n",((*(z+k))+1)*Ret,((*(z+k))+1),(utaut_sl[k]*utaut_sl[k])/(utauavg*utauavg),(utaut_nsl[k]*utaut_nsl[k])/(utauavg*utauavg));
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"shear-rate-slip-no-slip-tau-total-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"z<sup>+</sup>\", \"z/h\", \"z<sup>+ls</sup>\", \"z<sup>+lns</sup>\"");
        fprintf(tcp,", \"S<sub>sl</sub><sup>+l</sup>\", \"S<sub>ns</sub><sup>+l</sup>\"\n");
        for (k=1;k<(up-1);k++)
        {
                fprintf(tcp,"%f %f ",((*(z+k))+1)*Ret,((*(z+k))+1));
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret*utaut_sl[k]/utauavg);
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret*utaut_nsl[k]/utauavg);
		fprintf(tcp,"%.9f ",dudxsl[k]*nu/(utaut_sl[k]*utaut_sl[k]));
		fprintf(tcp,"%.9f\n",dudxnsl[k]*nu/(utaut_nsl[k]*utaut_nsl[k]));
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"uavg-test-tau-total-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	for (k=1;k<(up-1);k++)
		fprintf(tcp,"%f %f\n",((*(z+k))+1)*Ret*utaut_nsl[k]/utauavg,0.5*(U1[k]+U1[zd-1-k])/utaut_nsl[k]+integ[k]/(0.4*utaut_nsl[k]));
	fclose(tcp);

//	utns=dudxnsl[0]*nu;
//	printf("new ut_ns=%f\n",sqrt(utns));
//
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			if (SS[j] == 5)
			{
				uu1np[k] += (utin[j][k]/wstrs[j]);
                        	uv1np[k] += (uv[j][k]/wstrs[j]);
                        	uw1np[k] += (uw[j][k]/wstrs[j]);
                        	vv1np[k] += (vtin[j][k]/wstrs[j]);
                        	vw1np[k] += (vw[j][k]/wstrs[j]);
                        	ww1np[k] += (wtin[j][k]/wstrs[j]);
                        	UW1np[k] += (UW[j][k]/wstrs[j]);

                        	U1np[k] += (UA[j][k]/sqrt(wstrs[j]));

				dAVG1np[k] += (dAVG[j][k]/wstrs[j]);
                        	dfluc1np[k] += (dfluc[j][k]/wstrs[j]);
			}
		}
	}
	
	printf("check ct noslip:%d\n",ct);
	for (k=0;k<zd;k++)
	{
		uu1np[k] /= ((double)ct);
                uv1np[k] /= ((double)ct);
                uw1np[k] /= ((double)ct);
                vv1np[k] /= ((double)ct);
                vw1np[k] /= ((double)ct);
                ww1np[k] /= ((double)ct);
                UW1np[k] /= ((double)ct);
                U1np[k] /= ((double)ct);

                dAVG1np[k] /= ((double)ct);
                dfluc1np[k] /= ((double)ct);
	}

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

	sprintf(fn,"Aniso-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"II\",\"III\"\n");
	for (k=up-1;k>0;k--)
		fprintf(tcp,"%f %.9f %.9f\n",((*(z+k))+1)*Ret,IIa[k],IIIa[k]);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-no-slip-collected-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\",\"<greek>R</greek><sup>+</sup>\",\"<greek>r</greek><sup>+</sup>\"\n");
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
		fprintf(tcp,"%f ",0.5*(UW1[k]-UW1[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(dAVG1[k]+dAVG1[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",sqrt(0.5*(dfluc1[k]+dfluc1[zd-1-k]))/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-no-slip-collected-local-pl-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
        fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\",\"<greek>R</greek><sup>+</sup>\",\"<greek>r</greek><sup>+</sup>\"\n");
	for (k=0;k<up;k++)
        {
                fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
                fprintf(tcp,"%f ",((*(z+k))+1));
                fprintf(tcp,"%f ",0.5*(U1np[k]+U1np[zd-1-k]));
                fprintf(tcp,"%f ",sqrt(0.5*(uu1np[k]+uu1np[zd-1-k])));
                fprintf(tcp,"%f ",sqrt(0.5*(vv1np[k]+vv1np[zd-1-k])));
                fprintf(tcp,"%f ",sqrt(0.5*(ww1np[k]+ww1np[zd-1-k])));
                fprintf(tcp,"%f ",0.5*(uv1np[k]-uv1np[zd-1-k]));
                fprintf(tcp,"%f ",0.5*(uw1np[k]-uw1np[zd-1-k]));
                fprintf(tcp,"%f ",0.5*(vw1np[k]-vw1np[zd-1-k]));
                fprintf(tcp,"%f ",0.5*(UW1np[k]-UW1np[zd-1-k]));
                fprintf(tcp,"%f ",0.5*(dAVG1np[k]+dAVG1np[zd-1-k]));
                fprintf(tcp,"%f\n",sqrt(0.5*(dfluc1np[k]+dfluc1np[zd-1-k])));
        }
        fprintf(tcp,"\n");
        fclose(tcp);

// first rms then avg
	sprintf(fn,"stat-no-slip-collected-pl-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\",\"<greek>R</greek><sup>+</sup>\",\"<greek>r</greek><sup>+</sup>\"\n");
//	fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
	for (k=0;k<up;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",0.5*(U1[k]+U1[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",0.5*(uusq1[k]+uusq1[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",0.5*(vvsq1[k]+vvsq1[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",0.5*(wwsq1[k]+wwsq1[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",0.5*(uv1[k]-uv1[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(uw1[k]-uw1[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(vw1[k]-vw1[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(UW1[k]-UW1[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(dAVG1[k]+dAVG1[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",0.5*(dflucsq1[k]+dflucsq1[zd-1-k])/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	utns=sqrt(utns);

	sprintf(fn,"stat-no-slip-collected-w-its-ut-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
        fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"");
        fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
        for (k=0;k<up;k++)
        {
                fprintf(tcp,"%f ",((*(z+k))+1)*Ret*utns/utauavg);
                fprintf(tcp,"%f ",((*(z+k))+1));
                fprintf(tcp,"%f ",0.5*(U1[k]+U1[zd-1-k])/utns);
                fprintf(tcp,"%f ",sqrt(0.5*(uu1[k]+uu1[zd-1-k]))/utns);
                fprintf(tcp,"%f ",sqrt(0.5*(vv1[k]+vv1[zd-1-k]))/utns);
                fprintf(tcp,"%f ",sqrt(0.5*(ww1[k]+ww1[zd-1-k]))/utns);
                fprintf(tcp,"%f ",0.5*(uv1[k]-uv1[zd-1-k])/(utns*utns));
                fprintf(tcp,"%f ",0.5*(uw1[k]-uw1[zd-1-k])/(utns*utns));
                fprintf(tcp,"%f ",0.5*(vw1[k]-vw1[zd-1-k])/(utns*utns));
                fprintf(tcp,"%f ",0.5*(UW1[k]-UW1[zd-1-k])/(utns*utns));
		fprintf(tcp,"%f ",nu*(dudxnsl[k])/(utns*utns));
		fprintf(tcp,"%f\n",-0.5*(uw1[k]-uw1[zd-1-k])/(utns*utns)-0.5*(UW1[k]-UW1[zd-1-k])/(utns*utns)+nu*(dudxnsl[k])/(utns*utns));
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"stat-no-slip-collected-w-tau-total-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES =\"Z+\",\"Z+t\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
        fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"");
        fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
        for (k=1;k<(up-1);k++)
        {
                fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret*utaut_nsl[k]/utauavg);
                fprintf(tcp,"%f ",((*(z+k))+1));
                fprintf(tcp,"%f ",0.5*(U1[k]+U1[zd-1-k])/utaut_nsl[k]);
                fprintf(tcp,"%f ",sqrt(0.5*(uu1[k]+uu1[zd-1-k]))/utaut_nsl[k]);
                fprintf(tcp,"%f ",sqrt(0.5*(vv1[k]+vv1[zd-1-k]))/utaut_nsl[k]);
                fprintf(tcp,"%f ",sqrt(0.5*(ww1[k]+ww1[zd-1-k]))/utaut_nsl[k]);
                fprintf(tcp,"%f ",0.5*(uv1[k]-uv1[zd-1-k])/(utaut_nsl[k]*utaut_nsl[k]));
                fprintf(tcp,"%f ",0.5*(uw1[k]-uw1[zd-1-k])/(utaut_nsl[k]*utaut_nsl[k]));
                fprintf(tcp,"%f ",0.5*(vw1[k]-vw1[zd-1-k])/(utaut_nsl[k]*utaut_nsl[k]));
                fprintf(tcp,"%f ",0.5*(UW1[k]-UW1[zd-1-k])/(utaut_nsl[k]*utaut_nsl[k]));
                fprintf(tcp,"%f ",nu*(dudxnsl[k])/(utaut_nsl[k]*utaut_nsl[k]));
                fprintf(tcp,"%f\n",-0.5*(uw1[k]-uw1[zd-1-k])/(utaut_nsl[k]*utaut_nsl[k])-0.5*(UW1[k]-UW1[zd-1-k])/(utaut_nsl[k]*utaut_nsl[k])+nu*(dudxnsl[k])/(utaut_nsl[k]*utaut_nsl[k]));
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"stat-no-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\",\"<greek>R</greek><sup>+</sup>\",\"<greek>r</greek><sup>+</sup>\"\n");
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
		fprintf(tcp,"%f ",(UW1[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(dAVG1[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",sqrt(dfluc1[k])/(utauavg*utauavg));
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

	sprintf(fn,"Aniso-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"II\",\"III\"\n");
	for (k=up-1;k>0;k--)
		fprintf(tcp,"%f %.9f %.9f\n",((*(z+k))+1)*Ret,IIa[k],IIIa[k]);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-slip-collected-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<greek>R</greek><sup>+</sup>\",\"<greek>r</greek><sup>+</sup>\",\"U<sup>+</sup>-U<sub>s</sub><sup>+</sup>\"\n");
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
		fprintf(tcp,"%f ",0.5*(UW2[k]-UW2[zd-1-k])/(utauavg*utauavg));
//		fprintf(tcp,"%f ",nu*(dudxsl[k])/(utauavg*utauavg));
//                fprintf(tcp,"%f\n",-0.5*(uw2[k]-uw2[zd-1-k])/(utauavg*utauavg)-0.5*(UW2[k]-UW2[zd-1-k])/(utauavg*utauavg)+nu*(dudxsl[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(dAVG2[k]+dAVG2[zd-1-k])/(utauavg*utauavg));
                fprintf(tcp,"%f ",sqrt(0.5*(dfluc2[k]+dfluc2[zd-1-k]))/(utauavg*utauavg));
		fprintf(tcp,"%f\n",0.5*(U2[k]+U2[zd-1-k])/utauavg - 0.5*(uavg[1]+uavg[zd-1-1])/utauavg + ((*(z+1))+1)*Ret);
	}
	fprintf(tcp,"\n");
	fclose(tcp);

// Data of slip stripes Non-dimensionalized with ut of the no-slip surface
	sprintf(fn,"stat-slip-collected-w-its-ut-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
        fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"");
        fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
        for (k=0;k<up;k++)
        {
                fprintf(tcp,"%f ",((*(z+k))+1)*Ret*utns/utauavg);
                fprintf(tcp,"%f ",((*(z+k))+1));
                fprintf(tcp,"%f ",0.5*(U2[k]+U2[zd-1-k])/utns);
                fprintf(tcp,"%f ",sqrt(0.5*(uu2[k]+uu2[zd-1-k]))/utns);
                fprintf(tcp,"%f ",sqrt(0.5*(vv2[k]+vv2[zd-1-k]))/utns);
                fprintf(tcp,"%f ",sqrt(0.5*(ww2[k]+ww2[zd-1-k]))/utns);
                fprintf(tcp,"%f ",0.5*(uv2[k]-uv2[zd-1-k])/(utns*utns));
                fprintf(tcp,"%f ",0.5*(uw2[k]-uw2[zd-1-k])/(utns*utns));
                fprintf(tcp,"%f ",0.5*(vw2[k]-vw2[zd-1-k])/(utns*utns));
                fprintf(tcp,"%f ",0.5*(UW2[k]-UW2[zd-1-k])/(utns*utns));
                fprintf(tcp,"%f ",nu*(dudxsl[k])/(utns*utns));
                fprintf(tcp,"%f\n",-0.5*(uw2[k]-uw2[zd-1-k])/(utns*utns)-0.5*(UW2[k]-UW2[zd-1-k])/(utns*utns)+nu*(dudxsl[k])/(utns*utns));
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"stat-slip-collected-w-tau-total-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES =\"Z+\",\"Z+t\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
        fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"");
        fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
        for (k=1;k<(up-1);k++)
        {
                fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret*utaut_sl[k]/utauavg);
                fprintf(tcp,"%f ",((*(z+k))+1));
                fprintf(tcp,"%f ",0.5*(U2[k]+U2[zd-1-k])/utaut_sl[k]);
                fprintf(tcp,"%f ",sqrt(0.5*(uu2[k]+uu2[zd-1-k]))/utaut_sl[k]);
                fprintf(tcp,"%f ",sqrt(0.5*(vv2[k]+vv2[zd-1-k]))/utaut_sl[k]);
                fprintf(tcp,"%f ",sqrt(0.5*(ww2[k]+ww2[zd-1-k]))/utaut_sl[k]);
                fprintf(tcp,"%f ",0.5*(uv2[k]-uv2[zd-1-k])/(utaut_sl[k]*utaut_sl[k]));
                fprintf(tcp,"%f ",0.5*(uw2[k]-uw2[zd-1-k])/(utaut_sl[k]*utaut_sl[k]));
                fprintf(tcp,"%f ",0.5*(vw2[k]-vw2[zd-1-k])/(utaut_sl[k]*utaut_sl[k]));
                fprintf(tcp,"%f ",0.5*(UW2[k]-UW2[zd-1-k])/(utaut_sl[k]*utaut_sl[k]));
                fprintf(tcp,"%f ",nu*(dudxsl[k])/(utaut_sl[k]*utaut_sl[k]));
                fprintf(tcp,"%f\n",-0.5*(uw2[k]-uw2[zd-1-k])/(utaut_sl[k]*utaut_sl[k])-0.5*(UW2[k]-UW2[zd-1-k])/(utaut_sl[k]*utaut_sl[k])+nu*(dudxsl[k])/(utaut_sl[k]*utaut_sl[k]));
        }
        fprintf(tcp,"\n");
        fclose(tcp);

//first rms then avg
	sprintf(fn,"stat-slip-collected-pl-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<greek>R</greek><sup>+</sup>\",\"<greek>r</greek><sup>+</sup>\"\n");
//	fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
	for (k=0;k<up;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",0.5*(U2[k]+U2[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",0.5*(uusq2[k]+uu2[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",0.5*(vvsq2[k]+vv2[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",0.5*(wwsq2[k]+ww2[zd-1-k])/utauavg);
		fprintf(tcp,"%f ",0.5*(uv2[k]-uv2[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(uw2[k]-uw2[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(vw2[k]-vw2[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(UW2[k]-UW2[zd-1-k])/(utauavg*utauavg));
//		fprintf(tcp,"%f ",nu*(dudxsl[k])/(utauavg*utauavg));
//                fprintf(tcp,"%f\n",-0.5*(uw2[k]-uw2[zd-1-k])/(utauavg*utauavg)-0.5*(UW2[k]-UW2[zd-1-k])/(utauavg*utauavg)+nu*(dudxsl[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(dAVG2[k]+dAVG2[zd-1-k])/(utauavg*utauavg));
                fprintf(tcp,"%f\n",0.5*(dflucsq2[k]+dflucsq2[zd-1-k])/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-slip-collected-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/h\",\"U<sup>+</sup>\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"\n");
//	fprintf(tcp,",\"dUdz<sup>+</sup>\",\"Ttotal\"\n");
	fprintf(tcp,",\"<greek>R</greek><sup>+</sup>\",\"<greek>r</greek><sup>+</sup>\"\n");
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
		fprintf(tcp,"%f ",(UW2[k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(dAVG2[k])/(utauavg*utauavg));
                fprintf(tcp,"%f\n",sqrt(dfluc2[k])/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

//average of rms
	sprintf(fn,"utin-2DAVG-rms.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(uusq3[k]+uusq3[zd-1-k])/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"vtin-2DAVG-rms.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(vvsq3[k]+vvsq3[zd-1-k])/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"wtin-2DAVG-rms.dat");
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(wwsq3[k]+wwsq3[zd-1-k])/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"dtin-2DAVG-rms.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"p<sup>'+</sup><sub>rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(dflucsq3[k]+dflucsq3[zd-1-k])/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"uv-abs-slip-noslip-2DAVG.dat",SLPTXT);
        tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES = \"z<sup>+</sup>\", \"z<sup>+<greek>*</greek></sup>\", \"z/h\", \"<uv><sup>+</sup><sub>nsl</sub>\", \"<uv><sup>+</sup><sub>sl</sub>\"");
	fprintf(tcp,", \"<uv><sup>+<greek>*</greek></sup><sub>nsl</sub>\", \"<uv><sup>+<greek>*</greek></sup><sub>sl</sub>\"\n");
	for (k=0;k<up;k++)
	{
		fprintf(tcp,"%f %f %f",((*(z+k))+1)*Ret,((*(z+k))+1)*Ret*utns/utauavg,((*(z+k))+1));
		fprintf(tcp," %f %f",0.5*(absuvnsl[k]+absuvnsl[zd-1-k])/(utauavg*utauavg),0.5*(absuvsl[k]+absuvsl[zd-1-k])/(utauavg*utauavg));
		fprintf(tcp," %f %f\n",0.5*(absuvnsl[k]+absuvnsl[zd-1-k])/(utns*utns),0.5*(absuvsl[k]+absuvsl[zd-1-k])/(utns*utns));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"dtin-slip-collected-2DAVG.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"p<sup>+'</sup><sub>rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	fprintf(tcp,"%f %f\n", 0., 0.);
        for (k=1;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(dfluc2[k]+dfluc2[zd-1-k]))/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"dtin-no-slip-collected-2DAVG.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"p<sup>+'</sup><sub>rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	fprintf(tcp,"%f %f\n", 0., 0.);
        for (k=1;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(dfluc1[k]+dfluc1[zd-1-k]))/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"davg-slip-collected-2DAVG.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<P><sup>+</sup>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(dAVG2[k]+dAVG2[zd-1-k])/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

        sprintf(fn,"davg-no-slip-collected-2DAVG.dat");
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<P><sup>+</sup>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(dAVG1[k]+dAVG1[zd-1-k])/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	free(z);
}
