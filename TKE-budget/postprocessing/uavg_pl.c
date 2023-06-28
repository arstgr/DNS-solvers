#include "definitions.h"

void uAVG_pl(int xd, int yd, int zd, int st, int en)
{
	int i,j,k,a,c;
	FILE *sv;
	FILE *tcp,*wstr;
	char fn[20];
	int num=(en-st)+1;
	DP *z,dplus,Ret;
	DP utauavg=0.,utvt,time,nu;
	int cnt=0,ct=0,rd;
	int up,gst;
	DP ubavg=0.,ubt=0.,twall,val,tau,ub=0.;
	char stri[60];
	DP uttemp,uttemp1,ttemp;
	DP uAVG[yd][zd],uAVGt[yd][zd],vAVG[yd][zd],vAVGt[yd][zd],wAVG[yd][zd],wAVGt[yd][zd];
	DP oavg[yd][zd],oavgt[yd][zd];
	DP wstravg[yd],u1,u2,u3;
	DP uv2d[yd][zd],uv2dt[yd][zd],uw2d[yd][zd],uw2dt[yd][zd],uu2d[yd][zd],vv2d[yd][zd],ww2d[yd][zd],uu2dt[yd][zd],vv2dt[yd][zd],ww2dt[yd][zd];
	DP uv2dm[(gplus+wplus)][zd],uv2d2[yd][1+(zd-1)/2],uw2dm[(gplus+wplus)][zd],uw2d2[yd][1+(zd-1)/2];
	DP uu2dm[(gplus+wplus)][zd],uu2d2[yd][1+(zd-1)/2],vv2dm[(gplus+wplus)][zd],vv2d2[yd][1+(zd-1)/2],ww2dm[(gplus+wplus)][zd],ww2d2[yd][1+(zd-1)/2];
	DP UAVG[zd],VAVG[zd],WAVG[zd];

	DP Sd[yd],Su[yd],wstrr=0.;

	up=1+(zd-1)/2;
	num/=DT;
	z=(DP *)calloc(zd,sizeof(DP));
	gst=wst;

	for (j=0;j<yd;j++)
        {
                Sd[j] = -3.;
                Su[j] = 3.;
        }

        for (j=wst;j<yd;j+=(gplus+wplus))
        {
                for (a=0;a<gplus;a++)
                {
                        Sd[j+a] = -5.;
                        Su[j+a] = 5.;
                }
        }

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
	printf("uAVG: ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);

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
		UAVG[k] = 0.;
		VAVG[k] = 0.;
		WAVG[k] = 0.;
	}
	
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			uAVG[j][k]=0.;
			uAVGt[j][k]=0.;

			vAVG[j][k]=0.;
                        vAVGt[j][k]=0.;

			wAVG[j][k]=0.;
                        wAVGt[j][k]=0.;
			
			uv2d[j][k] = 0.;
			uw2d[j][k] = 0.;

			uu2d[j][k] = 0.;
			vv2d[j][k] = 0.;
			ww2d[j][k] = 0.;
		}
	}

	for (j=0;j<(gplus+wplus);j++)
	{
		for (k=0;k<zd;k++)
		{
			uv2dm[j][k] =0.;
			uw2dm[j][k] =0.;

			uu2dm[j][k] =0.;
			vv2dm[j][k] =0.;
			ww2dm[j][k] =0.;
		}
	}
	for (j=0;j<yd;j++)
	{
		for (k=0;k<up;k++)
		{
			uv2d2[j][k]=0.;
			uw2d2[j][k]=0.;

			uu2d2[j][k]=0.;
			vv2d2[j][k]=0.;
			ww2d2[j][k]=0.;
		}
	}
			
	
	for (i=st;i<(en+1);i+=DT)
	{
		sprintf(fn,"uAVG-2d.%.3d.%.4d",0,i);
		sv=fopen(fn,"rb");
		fread(&uAVGt[0][0],sizeof(DP),yd*zd,sv);
		fclose(sv);
		
		sprintf(fn,"vAVG-2d.%.3d.%.4d",0,i);
                sv=fopen(fn,"rb");
                fread(&vAVGt[0][0],sizeof(DP),yd*zd,sv);
                fclose(sv);

		sprintf(fn,"wAVG-2d.%.3d.%.4d",0,i);
                sv=fopen(fn,"rb");
                fread(&wAVGt[0][0],sizeof(DP),yd*zd,sv);
                fclose(sv);
		
		sprintf(fn,"vel-field-2d.%.3d.%.4d",0,i);
                sv=fopen(fn,"rb");
//		fseek(sv,3*yd*zd*sizeof(DP),SEEK_SET);
		fread(&uu2dt[0][0],sizeof(DP),yd*zd,sv);
		fread(&vv2dt[0][0],sizeof(DP),yd*zd,sv);
		fread(&ww2dt[0][0],sizeof(DP),yd*zd,sv);
                fread(&uv2dt[0][0],sizeof(DP),yd*zd,sv);
		fread(&uw2dt[0][0],sizeof(DP),yd*zd,sv);
                fclose(sv);
                
                sprintf(fn,"Mean-X-Vorticity-2d.%.3d.%.4d",0,i);
                sv=fopen(fn,"rb");
                fread(&oavgt[0][0],sizeof(DP),yd*zd,sv);
                fclose(sv);

		for (j=0;j<yd;j++)
		{
			for (k=0;k<zd;k++)
			{
				uAVG[j][k]+=uAVGt[j][k];
				vAVG[j][k]+=vAVGt[j][k];
				wAVG[j][k]+=wAVGt[j][k];
				oavg[j][k]+=oavgt[j][k];
				
				uv2d[j][k] += uv2dt[j][k];
				uw2d[j][k] += uw2dt[j][k];

				uu2d[j][k] += uu2dt[j][k];
				vv2d[j][k] += vv2dt[j][k];
				ww2d[j][k] += ww2dt[j][k];
			}
		}
	}
	
	for (j=0;j<yd;j++)
	{
		for (k=0;k<zd;k++)
		{
			uAVG[j][k]/=(num);
			vAVG[j][k]/=(num);
			wAVG[j][k]/=(num);
			oavg[j][k]/=(num);
			
			uv2d[j][k] /=(num);
			uw2d[j][k] /=(num);

			uu2d[j][k] /=(num);
			vv2d[j][k] /=(num);
			ww2d[j][k] /=(num);
		}
	}
	for (k=0;k<zd;k++)
	{
		for (j=0;j<yd;j++)
		{
			UAVG[k] += uAVG[j][k];
			VAVG[k] += vAVG[j][k];
			WAVG[k] += wAVG[j][k];
		}
	}
	for (k=0;k<zd;k++)
	{
		UAVG[k] /= yd;
		VAVG[k] /= yd;
		WAVG[k] /= yd;
	}

//*****************************************************************************************************
	for (j=gst;j<yd;j+=(gplus+wplus))
        {
                for (a=0;(a<(gplus+wplus))&&((a+j)<yd);a++)
                {
                        for (k=0;k<zd;k++)
                        {
                                c=(a+j-gst)%(gplus+wplus);
                                uv2dm[c][k]+=uv2d[a+j][k];
				uw2dm[c][k]+=uw2d[a+j][k];

				uu2dm[c][k]+=uu2d[a+j][k];
				vv2dm[c][k]+=vv2d[a+j][k];
				ww2dm[c][k]+=ww2d[a+j][k];
                        }
                }
        }

        for (a=0;a<gst;a++)
        {
                for (k=0;k<zd;k++)
                {
                        c=a+gst+gplus;
                        uv2dm[c][k]+=uv2d[a][k];
			uw2dm[c][k]+=uw2d[a][k];

			uu2dm[c][k]+=uu2d[a][k];
			vv2dm[c][k]+=vv2d[a][k];
			ww2dm[c][k]+=ww2d[a][k];
                }
        }

        c=yd/(gplus+wplus);
        for (a=0;a<(gplus+wplus);a++)
        {
                for (k=0;k<zd;k++)
                {
                        uv2dm[a][k]/=(double)c;
			uw2dm[a][k]/=(double)c;

			uu2dm[a][k]/=(double)c;
			vv2dm[a][k]/=(double)c;
			ww2dm[a][k]/=(double)c;
                }
        }

	for (j=gst;j<yd;j+=(gplus+wplus))
        {
                for (a=0;(a<(gplus+wplus))&&((a+j)<yd);a++)
                {
                        for (k=0;k<zd;k++)
                        {
                                c=(a+j-gst)%(gplus+wplus);
                                uv2d[j+a][k]=uv2dm[c][k];
				uw2d[j+a][k]=uw2dm[c][k];

				uu2d[j+a][k]=uu2dm[c][k];
				vv2d[j+a][k]=vv2dm[c][k];
				ww2d[j+a][k]=ww2dm[c][k];
                        }
                }
        }

        for (a=0;a<gst;a++)
        {
                for (k=0;k<zd;k++)
                {
                        c=a+gst+gplus;
                        uv2d[a][k]=uv2dm[c][k];
			uw2d[a][k]=uw2dm[c][k];

			uu2d[a][k]=uu2dm[c][k];
			vv2d[a][k]=vv2dm[c][k];
			ww2d[a][k]=ww2dm[c][k];
                }
        }

        for (j=0;j<yd;j++)
        {
                for (k=1;k<up;k++)
                {
                        uv2d2[j][k]=0.5*((uv2d[j][k])+(uv2d[j][zd-1-k]));
			uw2d2[j][k]=0.5*((uw2d[j][k])+(uw2d[j][zd-1-k]));

			uu2d2[j][k]=0.5*((uu2d[j][k])+(uu2d[j][zd-1-k]));
			vv2d2[j][k]=0.5*((vv2d[j][k])+(vv2d[j][zd-1-k]));
			ww2d2[j][k]=0.5*((ww2d[j][k])+(ww2d[j][zd-1-k]));
                }
        }

        for (j=0;j<yd;j++)
        {
                for (k=1;k<up;k++)
                {
                        uv2d[j][k]=uv2d2[j][k];
                        uv2d[j][zd-1-k]=uv2d2[j][k];
			
			uw2d[j][k]=uw2d2[j][k];
                        uw2d[j][zd-1-k]=uw2d2[j][k];

			uu2d[j][k]=uu2d2[j][k];
                        uu2d[j][zd-1-k]=uu2d2[j][k];

			vv2d[j][k]=vv2d2[j][k];
                        vv2d[j][zd-1-k]=vv2d2[j][k];

			ww2d[j][k]=ww2d2[j][k];
                        ww2d[j][zd-1-k]=ww2d2[j][k];
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
	wstrr = 0.;
	for (j=0;j<yd;j++)
                wstrr += wstravg[j];

	wstrr /= ((double)yd);

	sprintf(fn,"wstr-avg-unmod-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\", \"2y/(g+w)\",\"y/(<greek>n</greek>h/U<sub>b</sub>)<sup>1/2</sup>\", \"<greek>t</greek><sub>w</sub><sup>+0</sup>\"");
        fprintf(tcp,",\"<greek>t</greek>/<greek>r</greek>U<sub>b</sub><sup>2</sup>\",\"<<greek>t</greek>><sup><greek>*</greek></sup>/<greek>r</greek>U<sub>b</sub><sup>2</sup>\"");
        fprintf(tcp,",\"<<greek>t</greek>><sub>b</sub>/<greek>r</greek>U<sub>b</sub><sup>2</sup>\"");
        fprintf(tcp,",\"<<greek>t</greek>>/<greek>r</greek>U<sub>b</sub><sup>2</sup>\"");
        fprintf(tcp,",\"Su\",\"Sd\"\n");
        fprintf(tcp,"ZONE T=\"%s tau_w\", I=%d, F=POINT\n",stri,((*(z+1))+1)*Ret,yd);
        for (j=0;j<yd;j++)
        {
                fprintf(tcp,"%f %f %f %.9f %.9f",(j+0.5)*dplus,2.*(j+0.5)/(gplus+wplus),(j+0.5)/sqrt(nu*0.5*(zd-2.)/ub),wstravg[j]/(utbase*utbase),wstravg[j]/(ub*ub));
                fprintf(tcp," %.9f %.9f %.9f",(utauavg*utauavg)*((gplus+wplus)/(wplus))/(ub*ub),utbase*utbase/(ub*ub),(utauavg*utauavg)/(ub*ub));
                fprintf(tcp," %f %f\n",Su[j],Sd[j]);
        }
        fprintf(tcp,"\n");
        fclose(tcp);	

	for (j=wst;j<yd;j+=(gplus+wplus))
        {
                for (a=0;a<gplus;a++)
                {
                        wstravg[j+a] = 0.;
                }
        }
	wstrr = 0.;
	for (j=0;j<yd;j++)
		wstrr += wstravg[j];

	wstrr /= ((double)yd);
	for (j=0;j<yd;j++)
		wstravg[j] += (((utauavg*utauavg) - wstrr)*((gplus+wplus)/(wplus)));

	for (j=wst;j<yd;j+=(gplus+wplus))
        {
                for (a=0;a<gplus;a++)
                {
                        wstravg[j+a] = 0.;
                }
        }
	
	sprintf(stri,TEXT);

	sprintf(fn,"wstr-avg-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\", \"2y/(g+w)\",\"y/(<greek>n</greek>h/U<sub>b</sub>)<sup>1/2</sup>\", \"<greek>t</greek><sub>w</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<greek>t</greek>/<greek>r</greek>U<sub>b</sub><sup>2</sup>\",\"<<greek>t</greek>><sup><greek>*</greek></sup>/<greek>r</greek>U<sub>b</sub><sup>2</sup>\"");
	fprintf(tcp,",\"<<greek>t</greek>><sub>b</sub>/<greek>r</greek>U<sub>b</sub><sup>2</sup>\"");
	fprintf(tcp,",\"<<greek>t</greek>>/<greek>r</greek>U<sub>b</sub><sup>2</sup>\"");
	fprintf(tcp,",\"Su\",\"Sd\"\n");
        fprintf(tcp,"ZONE T=\"%s tau_w\", I=%d, F=POINT\n",stri,((*(z+1))+1)*Ret,yd);
        for (j=0;j<yd;j++)
	{
                fprintf(tcp,"%f %f %f %.9f %.9f",(j+0.5)*dplus,2.*(j+0.5)/(gplus+wplus),(j+0.5)/sqrt(nu*0.5*(zd-2.)/ub),wstravg[j]/(utauavg*utauavg),wstravg[j]/(ub*ub));
		fprintf(tcp," %.9f %.9f %.9f",(utauavg*utauavg)*((gplus+wplus)/(wplus))/(ub*ub),utbase*utbase/(ub*ub),(utauavg*utauavg)/(ub*ub));
		fprintf(tcp," %f %f\n",Su[j],Sd[j]);
	}
        fprintf(tcp,"\n");
        fclose(tcp);

//  Normalized in plus units *******************************************************************
	sprintf(fn,"uAVG-z1-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Y<sup>+</sup>\", \"U<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,((*(z+1))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",(j+0.5)*dplus,uAVG[j][1]/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"uAVG-z3-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Y<sup>+</sup>\", \"U<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,((*(z+2))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",(j+0.5)*dplus,uAVG[j][2]/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"uAVG-z5-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Y<sup>+</sup>\", \"U<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,((*(z+3))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",(j+0.5)*dplus,uAVG[j][3]/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"uAVG-z7-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Y<sup>+</sup>\", \"U<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,((*(z+4))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",(j+0.5)*dplus,uAVG[j][4]/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"uAVG-z15-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Y<sup>+</sup>\", \"U<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,((*(z+8))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",(j+0.5)*dplus,uAVG[j][8]/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"uAVG-z40-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Y<sup>+</sup>\", \"U<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,((*(z+21))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",(j+0.5)*dplus,uAVG[j][21]/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z60-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"Y<sup>+</sup>\", \"U<sup>+</sup>\"\n");
        fprintf(tcp,"ZONE T=\"%s z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,((*(z+41))+1)*Ret,yd);
        for (j=0;j<yd;j++)
                fprintf(tcp,"%f %.9f\n",(j+0.5)*dplus,uAVG[j][41]/utauavg);
        fprintf(tcp,"\n");
        fclose(tcp);

//  Normalized with bulk velocity and half width of the shear free ridge *****************************
	sprintf(fn,"uAVG-z1-gm-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(1-0.5)/((double)gplus),((*(z+1))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5)/((double)gplus),uAVG[j][1]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z3-gm-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(2-0.5)/((double)gplus),((*(z+2))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5)/((double)gplus),uAVG[j][2]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z5-gm-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(3-0.5)/((double)gplus),((*(z+3))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5)/((double)gplus),uAVG[j][3]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z7-gm-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(4-0.5)/((double)gplus),((*(z+4))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5)/((double)gplus),uAVG[j][4]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z15-gm-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(8-0.5)/((double)gplus),((*(z+8))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5)/((double)gplus),uAVG[j][8]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z40-gm-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(21-0.5)/((double)gplus),((*(z+21))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5)/((double)gplus),uAVG[j][21]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z60-gm-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(41-0.5)/((double)gplus),((*(z+41))+1)*Ret,yd);
        for (j=0;j<yd;j++)
                fprintf(tcp,"%f %.9f\n",2.*(j+0.5)/((double)gplus),uAVG[j][41]/ub);
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"uAVG-zw-gm-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(gplus-0.5)/((double)gplus),((*(z+gplus))+1)*Ret,yd);
	for (j=0;j<yd;j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5)/((double)gplus),uAVG[j][gplus]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

// Only one slice *****************************************************************
	sprintf(fn,"uAVG-z1-oner-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(1-0.5)/((double)gplus),((*(z+1))+1)*Ret,(gplus+wplus));
	for (j=wst;j<(wplus+gplus+wst);j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5-(gplus+wst))/((double)gplus),uAVG[j][1]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z3-oner-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(2-0.5)/((double)gplus),((*(z+2))+1)*Ret,(gplus+wplus));
	for (j=wst;j<(wplus+gplus+wst);j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5-(gplus+wst))/((double)gplus),uAVG[j][2]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z5-oner-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(3-0.5)/((double)gplus),((*(z+3))+1)*Ret,(gplus+wplus));
	for (j=wst;j<(wplus+gplus+wst);j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5-(gplus+wst))/((double)gplus),uAVG[j][3]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z7-oner-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(4-0.5)/((double)gplus),((*(z+4))+1)*Ret,(gplus+wplus));
	for (j=wst;j<(wplus+gplus+wst);j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5-(gplus+wst))/((double)gplus),uAVG[j][4]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z15-oner-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(8-0.5)/((double)gplus),((*(z+8))+1)*Ret,(gplus+wplus));
	for (j=wst;j<(wplus+gplus+wst);j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5-(gplus+wst))/((double)gplus),uAVG[j][8]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z40-oner-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(21-0.5)/((double)gplus),((*(z+21))+1)*Ret,(gplus+wplus));
	for (j=wst;j<(wplus+gplus+wst);j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5-(gplus+wst))/((double)gplus),uAVG[j][21]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uAVG-z60-oner-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(41-0.5)/((double)gplus),((*(z+41))+1)*Ret,(gplus+wplus));
        for (j=wst;j<(wplus+gplus+wst);j++)
                fprintf(tcp,"%f %.9f\n",2.*(j+0.5-(gplus+wst))/((double)gplus),uAVG[j][41]/ub);
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"uAVG-zw-oner-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"2Y/w\", \"U/U<sub>b</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s 2Z/w=%f z<sup>+</sup>=%f\", I=%d, F=POINT\n",stri,2.*(gplus-0.5)/((double)gplus),((*(z+gplus))+1)*Ret,(gplus+wplus));
	for (j=wst;j<(wplus+gplus+wst);j++)
		fprintf(tcp,"%f %.9f\n",2.*(j+0.5-(gplus+wst))/((double)gplus),uAVG[j][gplus]/ub);
	fprintf(tcp,"\n");
	fclose(tcp);

//      2D contour of mean cross flow over all of the ridges
	sprintf(fn,"cross-flow-plus-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE = \"Cross flow contour\"\n");
	fprintf(tcp,"VARIABLES = \"Y<sup>+</sup>\", \"Z<sup>+</sup>\", \"U<sup>+</sup><sub>y</sub>\", \"U<sup>+</sup><sub>z</sub>\",\"<greek>W</greek><sup>+</sup><sub>x</sub>\"\n");
	fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

	for (k=1;k<(zd-1);k++)
	{
		for (j=0;j<yd;j++)
		{
			fprintf(tcp,"%f %f %f %f %f\n",(j+0.5)*utauavg/nu,(k-0.5)*utauavg/nu,vAVG[j][k]/utauavg,wAVG[j][k]/utauavg,oavg[j][k]*nu/(utauavg*utauavg));
		}
	}
	fclose(tcp);

	sprintf(fn,"cross-flow-plus-oner-half-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"Cross flow contour only one ridge\"\n");
        fprintf(tcp,"VARIABLES = \"Y<sup>+</sup>\", \"Z<sup>+</sup>\", \"U<sup>+</sup><sub>y</sub>\", \"U<sup>+</sup><sub>z</sub>\",\"<greek>W</greek><sup>+</sup><sub>x</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",(gplus+wplus),up-1);

        for (k=1;k<up;k++)
        {
                for (j=wst;j<(wplus+gplus+wst);j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",(j+0.5-(gplus+wst))*utauavg/nu,(k-0.5)*utauavg/nu,vAVG[j][k]/utauavg,wAVG[j][k]/utauavg,oavg[j][k]*nu/(utauavg*utauavg));
                }
        }
        fclose(tcp);

	tcp=fopen("test.dat","w");
	for (j=0;j<yd;j++)
	{
		fprintf(tcp,"%.12f %.12f %.12f %.12f %.12f %.12f\n",uAVG[j][1],uAVG[j][zd-2],vAVG[j][1],vAVG[j][zd-2],wAVG[j][1],wAVG[j][zd-2]);
	}
	fclose(tcp);

	sprintf(fn,"cross-flow-bulk-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"Cross flow contour\"\n");
        fprintf(tcp,"VARIABLES = \"2Y/w\", \"2Z/w\", \"U<sub>y</sub>/U<sub>b</sub>\", \"U<sub>z</sub>/U<sub>b</sub>\",\"<greek>W</greek>w/2U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",2.*(j+0.5)/((double)gplus),2.*(k-0.5)/((double)gplus),vAVG[j][k]/ub,wAVG[j][k]/ub,0.5*gplus*oavg[j][k]/ub);
                }
        }
        fclose(tcp);

	sprintf(fn,"cross-flow-bulk-g-fudged-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"Cross flow contour subtracted with mean velocity\"\n");
        fprintf(tcp,"VARIABLES = \"2y/(g+w)\", \"z/g\", \"U<sub>y</sub>/U<sub>b</sub>\", \"U<sub>z</sub>/U<sub>b</sub>\",\"<greek>W</greek>w/2U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",2.*(j+0.5)/((double)(gplus+wplus)),(k-0.5)/((double)gplus),(vAVG[j][k]-VAVG[k])/ub,(wAVG[j][k]-WAVG[k])/ub,0.5*gplus*oavg[j][k]/ub);
                }
        }
        fclose(tcp);

	sprintf(fn,"cross-flow-bulk-gw-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"Cross flow contour\"\n");
        fprintf(tcp,"VARIABLES = \"2Y/(w+g)\", \"2Z/(w+g)\", \"U<sub>y</sub>/U<sub>b</sub>\", \"U<sub>z</sub>/U<sub>b</sub>\",\"<greek>W</greek>w/2U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",2.*(j+0.5)/((double)(gplus+wplus)),2.*(k-0.5)/((double)(gplus+wplus)),vAVG[j][k]/ub,wAVG[j][k]/ub,0.5*gplus*oavg[j][k]/ub);
                }
        }
        fclose(tcp);

	sprintf(fn,"cross-flow-bulk-g-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"Cross flow contour\"\n");
        fprintf(tcp,"VARIABLES = \"2y/(w+g)\", \"z/g\", \"U<sub>y</sub>/U<sub>b</sub>\", \"U<sub>z</sub>/U<sub>b</sub>\",\"<greek>W</greek>w/2U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",2.*(j+0.5)/((double)(gplus+wplus)),(k-0.5)/((double)(gplus)),vAVG[j][k]/ub,wAVG[j][k]/ub,0.5*gplus*oavg[j][k]/ub);
                }
        }
        fclose(tcp);

	sprintf(fn,"cross-flow-w-and-ut-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"Cross flow contour\"\n");
        fprintf(tcp,"VARIABLES = \"Y/g\", \"Z/g\", \"U<sub>y</sub>/U<sub>b</sub>\", \"U<sub>z</sub>/U<sub>b</sub>\",\"<greek>W</greek>w/2U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",(j+0.5)/((double)gplus),(k-0.5)/((double)gplus),vAVG[j][k]/utauavg,wAVG[j][k]/utauavg,nu*oavg[j][k]/(utauavg*utauavg));
                }
        }
        fclose(tcp);

	sprintf(fn,"cross-flow-w-and-ut-gpluswplus-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"Cross flow contour\"\n");
        fprintf(tcp,"VARIABLES = \"Y/g+w\", \"Z/g+w\", \"U<sub>y</sub>/U<sub>b</sub>\", \"U<sub>z</sub>/U<sub>b</sub>\",\"<greek>W</greek>w/2U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",2.*(j+0.5)/((double)(gplus+wplus)),2.*(k-0.5)/((double)(gplus+wplus)),vAVG[j][k]/utauavg,wAVG[j][k]/utauavg,nu*oavg[j][k]/(utauavg*utauavg));
                }
        }
        fclose(tcp);

	sprintf(fn,"cross-flow-w-and-ut-gplus-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"Cross flow contour\"\n");
        fprintf(tcp,"VARIABLES = \"Y/g+w\", \"Z/g\", \"U<sub>y</sub>/U<sub>b</sub>\", \"U<sub>z</sub>/U<sub>b</sub>\",\"<greek>W</greek>w/2U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",yd,(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=0;j<yd;j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",2.*(j+0.5)/((double)(gplus+wplus)),(k-0.5)/((double)(gplus)),vAVG[j][k]/utauavg,wAVG[j][k]/utauavg,nu*oavg[j][k]/(utauavg*utauavg));
                }
        }
        fclose(tcp);

	sprintf(fn,"cross-flow-bulk-oner-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"Cross flow contour only one ridge\"\n");
        fprintf(tcp,"VARIABLES = \"2Y/w\", \"2Z/w\", \"U<sub>y</sub>/U<sub>b</sub>\", \"U<sub>z</sub>/U<sub>b</sub>\",\"<greek>W</greek>w/2U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",(gplus+wplus),(zd-2));

        for (k=1;k<(zd-1);k++)
        {
                for (j=wst;j<(wplus+gplus+wst);j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",2.*(j+0.5-(gplus+wst))/((double)gplus),2.*(k-0.5)/((double)gplus),vAVG[j][k]/ub,wAVG[j][k]/ub,0.5*gplus*oavg[j][k]/ub);
                }
        }
        fclose(tcp);

	sprintf(fn,"cross-flow-bulk-oner-half-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"TITLE = \"Cross flow contour only one ridge\"\n");
        fprintf(tcp,"VARIABLES = \"2Y/w\", \"2Z/w\", \"U<sub>y</sub>/U<sub>b</sub>\", \"U<sub>z</sub>/U<sub>b</sub>\",\"<greek>W</greek>w/2U<sub>b</sub>\"\n");
        fprintf(tcp,"ZONE I=%d, J=%d, F=POINT\n",(gplus+wplus),up-1);

        for (k=1;k<up;k++)
        {
                for (j=wst;j<(wplus+gplus+wst);j++)
                {
                        fprintf(tcp,"%f %f %f %f %f\n",2.*(j+0.5-(gplus+wst))/((double)gplus),2.*(k-0.5)/((double)gplus),vAVG[j][k]/ub,wAVG[j][k]/ub,0.5*gplus*oavg[j][k]/ub);
                }
        }
        fclose(tcp);
	
	sprintf(fn,"uAVG-all-oner-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"2Y/w\", \"Z1\", \"Z2\", \"Z4\", \"Z8\", \"Z21\"");
        fprintf(tcp,", \"Z41\", \"Z64\"\n");
        for (j=wst;j<(wplus+gplus+wst);j++)
        {
                fprintf(tcp,"%f ",2.*(j+0.5-(gplus+wst))/((double)gplus));
                fprintf(tcp,"%.9f ",uAVG[j][1]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][2]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][4]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][8]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][21]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][41]/ub);
                fprintf(tcp,"%.9f\n",uAVG[j][64]/ub);
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"uAVG-all-c-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
	fprintf(tcp,"\"2y/(g+w)\",");
	fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\",\"S<sub>u</sub>\",\"S<sub>d</sub>\"\n",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
        for (j=0;j<yd;j++)
        {
                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
		fprintf(tcp,"%f ",2*(j+0.5)/(gplus+wplus));
                fprintf(tcp,"%.9f ",uAVG[j][1]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][2]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][4]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][8]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][21]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][41]/ub);
                fprintf(tcp,"%.9f ",uAVG[j][64]/ub);
		fprintf(tcp,"%f %f\n",Su[j],Sd[j]);
        }
        fprintf(tcp,"\n");
        fclose(tcp);
	
	sprintf(fn,"uv2d-pl-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
//        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\", \"Z=1\", \"Z=2\", \"Z=4\", \"Z=8\", \"Z=21\"");
//        fprintf(tcp,", \"Z=41\", \"Z=64\"\n");
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
                fprintf(tcp,"%.9f ",uv2d[j][1]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][2]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][4]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][8]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][21]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][41]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][64]/(utauavg*utauavg));
		fprintf(tcp,"%f %f ",Su[j],Sd[j]);
		fprintf(tcp,"%0.9f\n",utauavg/ub);
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"uv2d-pl-aux-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
	fprintf(tcp,"\"S<sub>u</sub>\",\"S<sub>d</sub>\"\n");
	for (j=0;j<yd;j++)
        {
		fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
		fprintf(tcp,"%f %f\n",Su[j],Sd[j]);
	}
	fprintf(tcp,"\n");
        fclose(tcp);
	
	sprintf(fn,"uw2d-pl-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
        fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
	fprintf(tcp,", \"Su\", \"Sd\"\n");
        for (j=0;j<yd;j++)
        {
                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
		fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
                fprintf(tcp,"%.9f ",uw2d[j][1]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][2]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][4]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][8]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][21]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][41]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][64]/(utauavg*utauavg));
		fprintf(tcp,"%f %f\n",Su[j],Sd[j]);
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"u-rms-2d-pl-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
        fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
        fprintf(tcp,", \"Su\", \"Sd\"\n");
        for (j=0;j<yd;j++)
        {
                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
                fprintf(tcp,"%.9f ",sqrt(uu2d[j][1]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(uu2d[j][2]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(uu2d[j][4]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(uu2d[j][8]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(uu2d[j][21]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(uu2d[j][41]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(uu2d[j][64]/(utauavg*utauavg)));
                fprintf(tcp,"%f %f\n",Su[j],Sd[j]);
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"v-rms-2d-pl-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
        fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
        fprintf(tcp,", \"Su\", \"Sd\"\n");
        for (j=0;j<yd;j++)
        {
                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
                fprintf(tcp,"%.9f ",sqrt(vv2d[j][1]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(vv2d[j][2]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(vv2d[j][4]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(vv2d[j][8]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(vv2d[j][21]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(vv2d[j][41]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(vv2d[j][64]/(utauavg*utauavg)));
                fprintf(tcp,"%f %f\n",Su[j],Sd[j]);
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	sprintf(fn,"w-rms-2d-pl-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"y<sup>+</sup>\",\"y<sup>+0</sup>\",\"y/h\",");
        fprintf(tcp,"\"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*0.5/(gplus+wplus),2.*1.5/(gplus+wplus),2.*3.5/(gplus+wplus),2.*7.5/(gplus+wplus),2.*20.5/(gplus+wplus));
        fprintf(tcp,", \"2z/(g+w)=%.3f\", \"2z/(g+w)=%.3f\"",2.*40.5/(gplus+wplus),2.*63.5/(gplus+wplus));
        fprintf(tcp,", \"Su\", \"Sd\"\n");
        for (j=0;j<yd;j++)
        {
                fprintf(tcp,"%f ",(j+0.5)*utauavg/nu);
                fprintf(tcp,"%f ",(j+0.5)*utbase/nu);
                fprintf(tcp,"%f ",2.*(j+0.5)/(zd-2.));
                fprintf(tcp,"%.9f ",sqrt(ww2d[j][1]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(ww2d[j][2]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(ww2d[j][4]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(ww2d[j][8]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(ww2d[j][21]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(ww2d[j][41]/(utauavg*utauavg)));
                fprintf(tcp,"%.9f ",sqrt(ww2d[j][64]/(utauavg*utauavg)));
                fprintf(tcp,"%f %f\n",Su[j],Sd[j]);
        }
        fprintf(tcp,"\n");
        fclose(tcp);
	
	sprintf(fn,"uv2d-all-oner-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"2Y/w\", \"Z1\", \"Z2\", \"Z4\", \"Z8\", \"Z21\"");
        fprintf(tcp,", \"Z41\", \"Z64\"\n");
        for (j=wst;j<(wplus+gplus+wst);j++)
        {
                fprintf(tcp,"%f ",2.*(j+0.5-(gplus+wst))/((double)gplus));
                fprintf(tcp,"%.9f ",uv2d[j][1]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][2]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][4]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][8]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][21]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uv2d[j][41]/(utauavg*utauavg));
                fprintf(tcp,"%.9f\n",uv2d[j][64]/(utauavg*utauavg));
        }
        fprintf(tcp,"\n");
        fclose(tcp);
	
	sprintf(fn,"uw2d-all-oner-%s.dat",SLPTXT);
        tcp=fopen(fn,"w");
        fprintf(tcp,"VARIABLES = \"2Y/w\", \"Z1\", \"Z2\", \"Z4\", \"Z8\", \"Z21\"");
        fprintf(tcp,", \"Z41\", \"Z64\"\n");
        for (j=wst;j<(wplus+gplus+wst);j++)
        {
                fprintf(tcp,"%f ",2.*(j+0.5-(gplus+wst))/((double)gplus));
                fprintf(tcp,"%.9f ",uw2d[j][1]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][2]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][4]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][8]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][21]/(utauavg*utauavg));
                fprintf(tcp,"%.9f ",uw2d[j][41]/(utauavg*utauavg));
                fprintf(tcp,"%.9f\n",uw2d[j][64]/(utauavg*utauavg));
        }
        fprintf(tcp,"\n");
        fclose(tcp);

	free(z);
}