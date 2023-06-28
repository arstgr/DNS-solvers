#include "definitions.h"

void vort_Line_results(int xd, int yd, int zd, int st, int en, int tinst)
{
	int i,j,k,a,c;
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
	double wx1[zd],wy1[zd],wz1[zd];
	double wx2[zd],wy2[zd],wz2[zd];
	double wx3[zd],wy3[zd],wz3[zd];
	double ox[yd][zd],oy[yd][zd],oz[yd][zd];
	double oxsq[yd][zd],oysq[yd][zd],ozsq[yd][zd];
	double oxt[yd][zd],oyt[yd][zd],ozt[yd][zd];
	double oxi[yd][zd],oyi[yd][zd],ozi[yd][zd];
	double ox1[zd],oy1[zd],oz1[zd];
	double ox2[zd],oy2[zd],oz2[zd];
	double oxsq1[zd],oysq1[zd],ozsq1[zd];
	double oxsq2[zd],oysq2[zd],ozsq2[zd];
	double oxsq3[zd],oysq3[zd],ozsq3[zd];
	int gst;
	double wstravg[yd],u1,u2,u3,wstrr;
	double uAVG[yd][zd],uAVGt[yd][zd];

	DP Sd[yd],Su[yd];
	DP ox2d[yd][zd],ox2dm[(gplus+wplus)][zd],ox2d2[yd][1+(zd-1)/2];
	DP oy2d[yd][zd],oy2dm[(gplus+wplus)][zd],oy2d2[yd][1+(zd-1)/2];
	DP oz2d[yd][zd],oz2dm[(gplus+wplus)][zd],oz2d2[yd][1+(zd-1)/2];

	DP oxnp[zd],oynp[zd],oznp[zd];
	double utaut_nsl[zd],utaut_sl[zd],utaut_avg[zd];
        DP uttotal,dsl[zd],dnsl[zd],dtot[zd];
        DP dudxsl[zd],dudxnsl[zd];
	DP Uns[zd],Us[zd];
        DP UW[yd][zd],UWt[yd][zd],uw[yd][zd],uwt[yd][zd],uw1[zd],uw2[zd],UW2[zd],UW1[zd];

	for (j=0;j<yd;j++)
	{
		Su[j]=3.;
		Sd[j]=-3.;
	}

	for (j=wst;j<yd;j+=(gplus+wplus))
        {
                for (a=0;a<gplus;a++)
                {
                        Sd[j+a] = -5.;
                        Su[j+a] = 5.;
                }
        }

	gst=wst;
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

	for (k=0;k<zd;k++)
	{	
		wx1[k] = 0.;
		wy1[k] = 0.;
		wz1[k] = 0.;

		wx2[k] = 0.;
		wy2[k] = 0.;
		wz2[k] = 0.;
		
		wx3[k] = 0.;
		wy3[k] = 0.;
		wz3[k] = 0.;

		ox1[k] = 0.;
		oy1[k] = 0.;
		oz1[k] = 0.;

		ox2[k] = 0.;
		oy2[k] = 0.;
		oz2[k] = 0.;

		oxsq1[k] = 0.;
		oysq1[k] = 0.;
		ozsq1[k] = 0.;

		oxsq2[k] = 0.;
		oysq2[k] = 0.;
		ozsq2[k] = 0.;

		oxsq3[k] = 0.;
		oysq3[k] = 0.;
		ozsq3[k] = 0.;

		oxnp[k] = 0.;
		oynp[k] = 0.;
		oznp[k] = 0.;

		Uns[k] = 0.;
                Us[k] = 0.;

                uw1[k] = 0.;
                uw2[k] = 0.;

                UW1[k] = 0.;
                UW2[k] = 0.;
	}

	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			ox[j][k] = 0.;
			oy[j][k] = 0.;
			oz[j][k] = 0.;

			oxsq[j][k] = 0.;
			oysq[j][k] = 0.;
			ozsq[j][k] = 0.;

			ox2d[j][k] = 0.;
                        oy2d[j][k] = 0.;
                        oz2d[j][k] = 0.;

			uAVG[j][k] = 0.;
			uAVGt[j][k] = 0.;

			uw[j][k] = 0.;
			uwt[j][k] = 0.;
		
			UW[j][k] = 0.;
			UWt[j][k] = 0.;
		}
	}

	for (j=0;j<(gplus+wplus);j++)
        {
                for (k=0;k<zd;k++)
                {
                        ox2dm[j][k] =0.;
                        oy2dm[j][k] =0.;
			oz2dm[j][k] =0.;
                }
        }
        for (j=0;j<yd;j++)
        {
                for (k=0;k<up;k++)
                {
                        ox2d2[j][k]=0.;
                        oy2d2[j][k]=0.;
			oz2d2[j][k] =0.;
                }
	}


	for (k=st;k<(en+1);k+=DT)
	{
		if (SLPF==0)
			sprintf(fn,"vort-fluc-2d-field.%.3d.%.4d",0,k);
		else if (SLPF==-1)
			sprintf(fn,"vort-fluc-2d-field.%.3d.%.4d",0,k);
		else if (SLPF==-2)
			sprintf(fn,"vort-fluc-2d-field.%.3d.%.4d",0,k);
		else if (SLPF==-3)
			sprintf(fn,"vort-fluc-2d-field.%.3d.%.4d",0,k);
		sv=fopen(fn,"rb");
		fread(&oxt[0][0],sizeof(DP),(yd*zd),sv);
		fread(&oyt[0][0],sizeof(DP),(yd*zd),sv);
		fread(&ozt[0][0],sizeof(DP),(yd*zd),sv);
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
				ox[j][a] += oxt[j][a];
				oy[j][a] += oyt[j][a];
				oz[j][a] += ozt[j][a];

				uAVG[j][a] += uAVGt[j][a];

				uw[j][a] += uwt[j][a];
                                UW[j][a] += UWt[j][a];
			}
		}
	}

	if (SLPF==0)
		sprintf(fn,"vort-fluc-2d-field.%.3d.%.4d",0,TINST);
	else if (SLPF==-1)
		sprintf(fn,"vort-fluc-2d-field.%.3d.%.4d",0,TINST);
	else if (SLPF==-2)
		sprintf(fn,"vort-fluc-2d-field.%.3d.%.4d",0,TINST);
	else if (SLPF==-3)
		sprintf(fn,"vort-fluc-2d-field.%.3d.%.4d",0,TINST);
	sv=fopen(fn,"rb");
	fread(&oxi[0][0],sizeof(DP),(yd*zd),sv);
	fread(&oyi[0][0],sizeof(DP),(yd*zd),sv);
	fread(&ozi[0][0],sizeof(DP),(yd*zd),sv);
	fclose(sv);

	for (j=0;j<yd;j++)
	{	
		for (k=0;k<zd;k++)
		{
			ox[j][k] /= ((double)num);
			oy[j][k] /= ((double)num);
			oz[j][k] /= ((double)num);

			uAVG[j][k] /= ((double)num);

			uw[j][k] /= ((double)num);
                        UW[j][k] /= ((double)num);
		}
	}
	for (j=0;j<yd;j++)
	{	
		for (k=0;k<zd;k++)
		{
			oxsq[j][k] = sqrt(ox[j][k]);
			oysq[j][k] = sqrt(oy[j][k]);
			ozsq[j][k] = sqrt(oz[j][k]);
		}
	}
	for (j=0;j<yd;j++)
	{	
		for (k=0;k<zd;k++)
		{
			oxsq3[k] += oxsq[j][k];
			oysq3[k] += oysq[j][k];
			ozsq3[k] += ozsq[j][k];
		}
	}
		
	for (k=0;k<zd;k++)
	{
		oxsq3[k] /= ((double)yd);
		oysq3[k] /= ((double)yd);
		ozsq3[k] /= ((double)yd);
	}
	
	ct=0;
	for (j=0;j<yd;j+=(gplus+wplus))
	{
		for (k=0;k<zd;k++)
		{	
			wx1[k] += ox[j][k];
			wy1[k] += oy[j][k];
			wz1[k] += oz[j][k];
			
			wx2[k] += ox[j+wst][k];
			wy2[k] += oy[j+wst][k];
			wz2[k] += oz[j+wst][k];
			
			wx3[k] += ox[j+wst+gplus/2][k];
			wy3[k] += oy[j+wst+gplus/2][k];
			wz3[k] += oz[j+wst+gplus/2][k];
		}
		ct++;
	}

	for (k=0;k<zd;k++)
	{	
		wx1[k] /= ((double)ct);
		wy1[k] /= ((double)ct);
		wz1[k] /= ((double)ct);
		
		wx2[k] /= ((double)ct);
		wy2[k] /= ((double)ct);
		wz2[k] /= ((double)ct);
			
		wx3[k] /= ((double)ct);
		wy3[k] /= ((double)ct);
		wz3[k] /= ((double)ct);
	}

	ct=0;
	for (j=wst;j<yd;j+=(gplus+wplus))
	{
		for (a=0;a<gplus;a++)
		{
			for (k=0;k<zd;k++)
			{
				ox2[k] += ox[j+a][k];
				oy2[k] += oy[j+a][k];
				oz2[k] += oz[j+a][k];

				oxsq2[k] += oxsq[j+a][k];
				oysq2[k] += oysq[j+a][k];
				ozsq2[k] += ozsq[j+a][k];

				Us[k] += uAVG[j+a][k];
                                uw2[k] += uw[j+a][k];
                                UW2[k] += UW[j+a][k];
			}
			ct++;
		}
	}
	for (k=0;k<zd;k++)
	{
		ox2[k] /= ((double)ct);
		oy2[k] /= ((double)ct);
		oz2[k] /= ((double)ct);

		oxsq2[k] /= ((double)ct);
		oysq2[k] /= ((double)ct);
		ozsq2[k] /= ((double)ct);

		Us[k] /= ((double)ct);
                uw2[k] /= ((double)ct);
                UW2[k] /= ((double)ct);
	}
	printf("ct slip: %d\n",ct);
	ct=0;
	for (j=(wst+gplus);j<(yd-wst);j+=(gplus+wplus))
	{
		for (a=0;a<wplus;a++)
		{
			for (k=0;k<zd;k++)
			{
				ox1[k] += ox[j+a][k];
				oy1[k] += oy[j+a][k];
				oz1[k] += oz[j+a][k];

				oxsq1[k] += oxsq[j+a][k];
				oysq1[k] += oysq[j+a][k];
				ozsq1[k] += ozsq[j+a][k];

				Uns[k] += uAVG[j+a][k];
                                uw1[k] += uw[j+a][k];
                                UW1[k] += UW[j+a][k];
			}
			ct++;
		}
	}
	for (a=0;a<wst;a++)
	{
		for (k=0;k<zd;k++)
		{
			ox1[k] += ox[a][k];
			oy1[k] += oy[a][k];
			oz1[k] += oz[a][k];

			Uns[k] += uAVG[a][k];
                        uw1[k] += uw[a][k];
                        UW1[k] += UW[a][k];

			ox1[k] += ox[yd-wst+a][k];
			oy1[k] += oy[yd-wst+a][k];
			oz1[k] += oz[yd-wst+a][k];

			oxsq1[k] += oxsq[a][k];
			oysq1[k] += oysq[a][k];
			ozsq1[k] += ozsq[a][k];

			oxsq1[k] += oxsq[yd-wst+a][k];
			oysq1[k] += oysq[yd-wst+a][k];
			ozsq1[k] += ozsq[yd-wst+a][k];

			Uns[k] += uAVG[yd-wst+a][k];
                        uw1[k] += uw[yd-wst+a][k];
                        UW1[k] += UW[yd-wst+a][k];
		}
		ct +=2;
	}

	for (k=0;k<zd;k++)
	{
		ox1[k] /= ((double)ct);
		oy1[k] /= ((double)ct);
		oz1[k] /= ((double)ct);

		oxsq1[k] /= ((double)ct);
		oysq1[k] /= ((double)ct);
		ozsq1[k] /= ((double)ct);

		Uns[k] /= ((double)ct);
		uw1[k] /= ((double)ct);
		UW1[k] /= ((double)ct);
	}
	printf("ct no-slip: %d\n",ct);
//*************************************************************************************************************
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			ox2d[j][k] = ox[j][k];
			oy2d[j][k] = oy[j][k];
			oz2d[j][k] = oz[j][k];
		}
	}	

	for (j=gst;j<yd;j+=(gplus+wplus))
        {
                for (a=0;(a<(gplus+wplus))&&((a+j)<yd);a++)
                {
                        for (k=0;k<zd;k++)
                        {
                                c=(a+j-gst)%(gplus+wplus);
                                ox2dm[c][k]+=ox2d[a+j][k];
                                oy2dm[c][k]+=oy2d[a+j][k];
				oz2dm[c][k]+=oz2d[a+j][k];
                        }
                }
        }

        for (a=0;a<gst;a++)
        {
                for (k=0;k<zd;k++)
                {
                        c=a+gst+gplus;
                        ox2dm[c][k]+=ox2d[a][k];
                        oy2dm[c][k]+=oy2d[a][k];
			oz2dm[c][k]+=oz2d[a][k];
                }
        }

        c=yd/(gplus+wplus);
        for (a=0;a<(gplus+wplus);a++)
        {
                for (k=0;k<zd;k++)
                {
                        ox2dm[a][k]/=(double)c;
                        oy2dm[a][k]/=(double)c;
			oz2dm[a][k]/=(double)c;
                }
        }

        for (j=gst;j<yd;j+=(gplus+wplus))
        {
                for (a=0;(a<(gplus+wplus))&&((a+j)<yd);a++)
                {
                        for (k=0;k<zd;k++)
                        {
                                c=(a+j-gst)%(gplus+wplus);
                                ox2d[j+a][k]=ox2dm[c][k];
                                oy2d[j+a][k]=oy2dm[c][k];
				oz2d[j+a][k]=oz2dm[c][k];
                        }
                }
        }
	
	for (a=0;a<gst;a++)
        {
                for (k=0;k<zd;k++)
                {
                        c=a+gst+gplus;
                        ox2d[a][k]=ox2dm[c][k];
                        oy2d[a][k]=oy2dm[c][k];
			oz2d[a][k]=oz2dm[c][k];
                }
        }

        for (j=0;j<yd;j++)
        {
                for (k=1;k<up;k++)
                {
                        ox2d2[j][k]=0.5*((ox2d[j][k])+(ox2d[j][zd-1-k]));
                        oy2d2[j][k]=0.5*((oy2d[j][k])+(oy2d[j][zd-1-k]));
			oz2d2[j][k]=0.5*((oz2d[j][k])+(oz2d[j][zd-1-k]));
                }
        }

        for (j=0;j<yd;j++)
        {
                for (k=1;k<up;k++)
                {
                        ox2d[j][k]=ox2d2[j][k];
                        ox2d[j][zd-1-k]=ox2d2[j][k];

                        oy2d[j][k]=oy2d2[j][k];
                        oy2d[j][zd-1-k]=oy2d2[j][k];

			oz2d[j][k]=oz2d2[j][k];
                        oz2d[j][zd-1-k]=oz2d2[j][k];
                }
        }

//*****************************************************************************************************
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
/*      for (j=0;j<yd;j++)
 *      	wstrr += wstravg[j];
 *
 *      wstrr /= ((double)yd);
 *      for (j=0;j<yd;j++)
 * 	     wstravg[j] += (((utauavg*utauavg) - wstrr)*((gplus+wplus)/(wplus)));
 **/
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
                utaut_sl[k] = sqrt((nu*dudxsl[k] - 0.5*(uw2[k]-uw2[zd-1-k]))/(1.-((*(z+k))+1)));
                utaut_avg[k] = sqrt((gplus*utaut_sl[k]*utaut_sl[k]+wplus*utaut_nsl[k]*utaut_nsl[k])/(gplus+wplus));
//              printf("ut_tot_ns=%f  ut_tot_sl=%f ut_tot_t=%f\n",utaut_nsl[k],utaut_sl[k],utaut_avg[k]);
        }
	printf("utau=%f\n",utauavg);

        uttotal = utaut_avg[up-2];
        printf("uttotal=%f error=%f\n",uttotal,fabs((utauavg-uttotal)/utauavg)*100.);

/*	
	for (k=0;k<(up-1);k++)
	{
                dsl[k] = (uttotal*uttotal - utaut_sl[k]*utaut_sl[k])*gplus;
                dnsl[k] = (utaut_nsl[k]*utaut_nsl[k] - uttotal*uttotal)*wplus;
                dtot[k] = 0.5*(dsl[k]+dnsl[k]);
        }
        for (k=0;k<(up-1);k++)
        {
                utaut_nsl[k] = sqrt(uttotal*uttotal + (dtot[k]/((double)wplus)));
                utaut_sl[k] = sqrt(uttotal*uttotal - (dtot[k]/((double)gplus)));
	}
*/
//*******************************************************************************************************
	
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			if (Su[j] == 3)
			{
				oxnp[k] += (nu*sqrt(ox2d[j][k])/wstravg[j]);
				oynp[k] += (nu*sqrt(oy2d[j][k])/wstravg[j]);
				oznp[k] += (nu*sqrt(oz2d[j][k])/wstravg[j]);
			}
		}
	}
	for (k=0;k<zd;k++)
	{
		oxnp[k] /= ((double)ct);
		oynp[k] /= ((double)ct);
		oznp[k] /= ((double)ct);
	}

	sprintf(stri,TEXT);
//**********************************************************************************************************
	sprintf(fn,"omega-rms-x-pl-no-slip-collected-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,oxnp[k]);
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-y-pl-no-slip-collected-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,oynp[k]);
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-z-pl-no-slip-collected-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,oznp[k]);
        fprintf(tcp,"\n");
        fclose(tcp);
//*************************************************************************************************************	

	sprintf(fn,"omega-rms-x-pl-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(wx1[k]+wx1[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-x-h-no-slip-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(wx1[k]+wx1[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(wx1[k]+wx1[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-y-pl-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(wy1[k]+wy1[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-y-h-no-slip-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(wy1[k]+wy1[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(wy1[k]+wy1[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-z-pl-no-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(wz1[k]+wz1[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-z-h-no-slip-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(wz1[k]+wz1[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(wz1[k]+wz1[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);
//**************************************************************************************************************************
	sprintf(fn,"omega-rms-x-pl-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(ox1[k]+ox1[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
//rms: first rms then mean
	sprintf(fn,"omega-rms-x-pl-no-slip-collected-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(oxsq1[k]+oxsq1[zd-1-k])*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-x-h-no-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(ox1[k]+ox1[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(ox1[k]+ox1[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-y-pl-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(oy1[k]+oy1[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

//first rms then avg
	sprintf(fn,"omega-rms-y-pl-no-slip-collected-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(oysq1[k]+oysq1[zd-1-k])*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-y-h-no-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(oy1[k]+oy1[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(oy1[k]+oy1[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-z-pl-no-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(oz1[k]+oz1[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

//first rms then avg
	sprintf(fn,"omega-rms-z-pl-no-slip-collected-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(ozsq1[k]+ozsq1[zd-1-k])*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-z-h-no-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(oz1[k]+oz1[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(oz1[k]+oz1[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

//**************************************************************************************************************************
	sprintf(fn,"omega-rms-x-pl-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(wx2[k]+wx2[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-x-h-edge-slip-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(wx2[k]+wx2[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(wx2[k]+wx2[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-y-pl-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(wy2[k]+wy2[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-y-h-edge-slip-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(wy2[k]+wy2[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(wy2[k]+wy2[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-z-pl-edge-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(wz2[k]+wz2[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-z-h-edge-slip-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(wz2[k]+wz2[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(wz2[k]+wz2[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);
//********************************************************************************************************************
	sprintf(fn,"omega-rms-x-pl-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(wx3[k]+wx3[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-x-h-slip-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(wx3[k]+wx3[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(wx3[k]+wx3[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-y-pl-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(wy3[k]+wy3[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-y-h-slip-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(wy3[k]+wy3[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(wy3[k]+wy3[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-z-pl-slip-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(wz3[k]+wz3[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-z-h-slip-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(wz3[k]+wz3[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(wz3[k]+wz3[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);
//**************************************************************************************************************
	sprintf(fn,"omega-rms-x-pl-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(ox2[k]+ox2[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

// first rms then avg
	sprintf(fn,"omega-rms-x-pl-slip-collected-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(oxsq2[k]+oxsq2[zd-1-k])*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-x-h-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(ox2[k]+ox2[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(ox2[k]+ox2[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-y-pl-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(oy2[k]+oy2[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

// first rms then avg
	sprintf(fn,"omega-rms-y-pl-slip-collected-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(oysq2[k]+oysq2[zd-1-k])*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-y-h-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(oy2[k]+oy2[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(oy2[k]+oy2[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-z-pl-slip-collected-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*(oz2[k]+oz2[zd-1-k]))*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

// first rms then avg
	sprintf(fn,"omega-rms-z-pl-slip-collected-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(ozsq2[k]+ozsq2[zd-1-k])*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-z-h-slip-collected-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
        for (k=0;k<up;k++)
                fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*(oz2[k]+oz2[zd-1-k]))*nu/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
                fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*(oz2[k]+oz2[zd-1-k]))*nu/(utauavg*utauavg));
        fprintf(tcp,"\n");
        fclose(tcp);

//**************************************************************************************************************
	sprintf(fn,"contour-omega-rms-plus-plus-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"omega-rms contour\"\n");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\", \"z<sup>+</sup>\", \"<greek>w</greek><sup>+</sup><sub>x</sub>\", \"<greek>w</greek><sup>+</sup><sub>y</sub>\",\"<greek>w</greek><sup>+</sup><sub>z</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",(j+0.5)*utauavg/nu,(k-0.5)*utauavg/nu,sqrt(ox[j][k])*nu/(utauavg*utauavg),sqrt(oy[j][k])*nu/(utauavg*utauavg),sqrt(oz[j][k])*nu/(utauavg*utauavg));
                }
        }
        fclose(tcp);

	sprintf(fn,"contour-omega-rms-bulk-plus-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"omega-rms contour\"\n");
        fprintf(tcp,"VARIABLES = \"y/g\", \"z/g\", \"<greek>w</greek><sup>+</sup><sub>x</sub>\", \"<greek>w</greek><sup>+</sup><sub>y</sub>\",\"<greek>w</greek><sup>+</sup><sub>z</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",(j+0.5)/gplus,(k-0.5)/gplus,sqrt(ox[j][k])*nu/(utauavg*utauavg),sqrt(oy[j][k])*nu/(utauavg*utauavg),sqrt(oz[j][k])*nu/(utauavg*utauavg));
                }
        }
        fclose(tcp);

	sprintf(fn,"contour-omega-rms-bulk-bulk-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"omega-rms contour\"\n");
        fprintf(tcp,"VARIABLES = \"y/g\", \"z/g\", \"<greek>w</greek><sub>x, b</sub>\", \"<greek>w</greek><sub>y, b</sub>\",\"<greek>w</greek><sub>z, b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",(j+0.5)/gplus,(k-0.5)/gplus,sqrt(ox[j][k])*0.5*(zd-2)/ub,sqrt(oy[j][k])*0.5*(zd-2)/ub,sqrt(oz[j][k])*0.5*(zd-2)/ub);
                }
        }
        fclose(tcp);
//******************************************************************************************************************************
	sprintf(fn,"contour-omega-inst-plus-plus-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"omega-rms contour TINS=%d\"\n",TINST);
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\", \"z<sup>+</sup>\", \"<greek>w</greek><sup>+</sup><sub>x</sub>\", \"<greek>w</greek><sup>+</sup><sub>y</sub>\",\"<greek>w</greek><sup>+</sup><sub>z</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",(j+0.5)*utauavg/nu,(k-0.5)*utauavg/nu,sqrt(oxi[j][k])*nu/(utauavg*utauavg),sqrt(oyi[j][k])*nu/(utauavg*utauavg),sqrt(ozi[j][k])*nu/(utauavg*utauavg));
                }
        }
        fclose(tcp);

	sprintf(fn,"contour-omega-inst-bulk-plus-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"omega-rms contour TINS=%d\"\n",TINST);
        fprintf(tcp,"VARIABLES = \"y/g\", \"z/g\", \"<greek>w</greek><sup>+</sup><sub>x</sub>\", \"<greek>w</greek><sup>+</sup><sub>y</sub>\",\"<greek>w</greek><sup>+</sup><sub>z</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",(j+0.5)/gplus,(k-0.5)/gplus,sqrt(oxi[j][k])*nu/(utauavg*utauavg),sqrt(oyi[j][k])*nu/(utauavg*utauavg),sqrt(ozi[j][k])*nu/(utauavg*utauavg));
                }
        }
        fclose(tcp);

	sprintf(fn,"contour-omega-inst-bulk-bulk-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"omega-rms contour TINS=%d\"\n",TINST);
        fprintf(tcp,"VARIABLES = \"y/g\", \"z/g\", \"<greek>w</greek><sub>x, b</sub>\", \"<greek>w</greek><sub>y, b</sub>\",\"<greek>w</greek><sub>z, b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",(j+0.5)/gplus,(k-0.5)/gplus,sqrt(oxi[j][k])*0.5*(zd-2)/ub,sqrt(oyi[j][k])*0.5*(zd-2)/ub,sqrt(ozi[j][k])*0.5*(zd-2)/ub);
                }
        }
        fclose(tcp);

//********************************************************************************************************************
//1D first rms then avg
	sprintf(fn,"omega-rms-x-pl-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(oxsq3[k]+oxsq3[zd-1-k])*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-y-pl-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(oysq3[k]+oysq3[zd-1-k])*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"omega-rms-z-pl-rms-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>W</greek><sup>+</sup><sub>i, rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(ozsq3[k]+ozsq3[zd-1-k])*nu/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
//********************************************************************************************************************
	sprintf(fn,"omega-x-rms-2d-pl-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
        fprintf(tcp,"\"2y/(g+w)\",");
        fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\",\"S<sub>u</sub>\",\"S<sub>d</sub>\",\"u<sub>t</sub>/U<sub>b</sub>\"\n",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
        for (j=0;j<yd;j++)
        {
                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
                fprintf(tcp,"%f ",2*(j+0.5)/(gplus+wplus));
                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][1])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][2])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][4])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][8])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][21])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][41])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][64])/(utauavg*utauavg));
                fprintf(tcp,"%f %f ",Su[j],Sd[j]);
                fprintf(tcp,"%0.9f\n",utauavg/ub);
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-x-rms-2d-pl-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
        fprintf(tcp,"\"2y/(g+w)\",");
        fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\",\"S<sub>u</sub>\",\"S<sub>d</sub>\",\"u<sub>t</sub>/U<sub>b</sub>\"\n",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
        for (j=0;j<yd;j++)
        {
		if (Sd[j]==-3)
		{
	                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
	                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
	                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
	                fprintf(tcp,"%f ",2*(j+0.5)/(gplus+wplus));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][1])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][2])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][4])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][8])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][21])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][41])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][64])/(wstravg[j]));
	                fprintf(tcp,"%f %f ",Su[j],Sd[j]);
	                fprintf(tcp,"%0.9f\n",utauavg/ub);
		}
		else
		{
                        fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
	                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
	                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
	                fprintf(tcp,"%f ",2*(j+0.5)/(gplus+wplus));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][1])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][2])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][4])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][8])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][21])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][41])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(ox2d[j][64])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%f %f ",Su[j],Sd[j]);
	                fprintf(tcp,"%0.9f\n",utauavg/ub);
		}
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-y-rms-2d-pl-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
        fprintf(tcp,"\"2y/(g+w)\",");
        fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\",\"S<sub>u</sub>\",\"S<sub>d</sub>\",\"u<sub>t</sub>/U<sub>b</sub>\"\n",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
        for (j=0;j<yd;j++)
        {
                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
                fprintf(tcp,"%f ",2*(j+0.5)/(gplus+wplus));
                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][1])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][2])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][4])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][8])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][21])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][41])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][64])/(utauavg*utauavg));
                fprintf(tcp,"%f %f ",Su[j],Sd[j]);
                fprintf(tcp,"%0.9f\n",utauavg/ub);
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-y-rms-2d-pl-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
        fprintf(tcp,"\"2y/(g+w)\",");
        fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\",\"S<sub>u</sub>\",\"S<sub>d</sub>\",\"u<sub>t</sub>/U<sub>b</sub>\"\n",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
        for (j=0;j<yd;j++)
        {
		if (Sd[j]==-3)
		{
	                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
	                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
	                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
	                fprintf(tcp,"%f ",2*(j+0.5)/(gplus+wplus));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][1])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][2])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][4])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][8])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][21])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][41])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][64])/(wstravg[j]));
	                fprintf(tcp,"%f %f ",Su[j],Sd[j]);
	                fprintf(tcp,"%0.9f\n",utauavg/ub);
		}
		else
		{
                        fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
	                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
	                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
	                fprintf(tcp,"%f ",2*(j+0.5)/(gplus+wplus));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][1])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][2])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][4])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][8])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][21])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][41])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oy2d[j][64])/(utauavg*utauavg));
	                fprintf(tcp,"%f %f ",Su[j],Sd[j]);
	                fprintf(tcp,"%0.9f\n",utauavg/ub);
		}
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-z-rms-2d-pl-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
        fprintf(tcp,"\"2y/(g+w)\",");
        fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\",\"S<sub>u</sub>\",\"S<sub>d</sub>\",\"u<sub>t</sub>/U<sub>b</sub>\"\n",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
        for (j=0;j<yd;j++)
        {
                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
                fprintf(tcp,"%f ",2*(j+0.5)/(gplus+wplus));
                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][1])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][2])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][4])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][8])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][21])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][41])/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][64])/(utauavg*utauavg));
                fprintf(tcp,"%f %f ",Su[j],Sd[j]);
                fprintf(tcp,"%0.9f\n",utauavg/ub);
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-z-rms-2d-pl-local-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
        fprintf(tcp,"\"2y/(g+w)\",");
        fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\",\"S<sub>u</sub>\",\"S<sub>d</sub>\",\"u<sub>t</sub>/U<sub>b</sub>\"\n",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
        for (j=0;j<yd;j++)
        {
		if (Sd[j]==-3)
		{
	                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
	                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
	                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
	                fprintf(tcp,"%f ",2*(j+0.5)/(gplus+wplus));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][1])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][2])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][4])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][8])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][21])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][41])/(wstravg[j]));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][64])/(wstravg[j]));
	                fprintf(tcp,"%f %f ",Su[j],Sd[j]);
	                fprintf(tcp,"%0.9f\n",utauavg/ub);
		}
		else
		{
                        fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
	                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
	                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
	                fprintf(tcp,"%f ",2*(j+0.5)/(gplus+wplus));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][1])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][2])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][4])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][8])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][21])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][41])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%.9f ",nu*sqrt(oz2d[j][64])*(gplus/((double)(gplus+wplus)))/(utauavg*utauavg));
	                fprintf(tcp,"%f %f ",Su[j],Sd[j]);
	                fprintf(tcp,"%0.9f\n",utauavg/ub);
		}
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"omega-rms-slip-no-slip-collected-w-tau-total-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"Z+\", \"Z+t\",\"Z/h\"");
	fprintf(tcp,", \"Oxns\",\"Oyns\",\"Ozns\",\"Oxs\",\"Oys\",\"Ozs\"\n");
        for (k=1;k<(up-1);k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret*utaut_sl[k]/utauavg);
                fprintf(tcp,"%f ",((*(z+k))+1)*Ret*utaut_nsl[k]/utauavg);
                fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",sqrt(0.5*(ox1[k]+ox1[zd-1-k]))*nu/(utaut_nsl[k]*utaut_nsl[k]));
		fprintf(tcp,"%f ",sqrt(0.5*(oy1[k]+oy1[zd-1-k]))*nu/(utaut_nsl[k]*utaut_nsl[k]));
		fprintf(tcp,"%f ",sqrt(0.5*(oz1[k]+oz1[zd-1-k]))*nu/(utaut_nsl[k]*utaut_nsl[k]));
		fprintf(tcp,"%f ",sqrt(0.5*(ox2[k]+ox2[zd-1-k]))*nu/(utaut_sl[k]*utaut_sl[k]));
		fprintf(tcp,"%f ",sqrt(0.5*(oy2[k]+oy2[zd-1-k]))*nu/(utaut_sl[k]*utaut_sl[k]));
		fprintf(tcp,"%f \n",sqrt(0.5*(oz2[k]+oz2[zd-1-k]))*nu/(utaut_sl[k]*utaut_sl[k]));
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	free(z);
}
