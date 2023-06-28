#include "definitions.h"

void time_stat(int xd, int yd, int zd, int st, int en)
{
	int i,j,k,a;
	FILE *sv;
	FILE *tcp,*wstr;
	DP *uavg,*vavg,*wavg,*davg;
	DP *utin,*vtin,*wtin,*reys,*dtin;
	DP *uavgt,*vavgt,*wavgt,*davgt;
	DP *utint,*vtint,*wtint,*reyst,*dtint;
	DP *sku,*skv,*skw,*skut,*skvt,*skwt;
	DP *skuw,*skp,*skuwt,*skpt;
	DP *kuu,*kuv,*kuw,*kuut,*kuvt,*kuwt;
	DP *kuuw,*kup,*kuuwt,*kupt;
	DP *kuuv,*kuuvt,*kuvw,*kuvwt;
	DP *skuv,*skuvt,*skvw,*skvwt;
	DP *uv,*uvt,*vw,*vwt;
	DP *UW,*UWt;
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
	DP utitsown;

	up=1+(zd-1)/2;
	num/=DT;
	z=(DP *)calloc(zd,sizeof(DP));

	uavg = (DP *)calloc(zd,sizeof(DP));
	vavg = (DP *)calloc(zd,sizeof(DP));
	wavg = (DP *)calloc(zd,sizeof(DP));
	davg = (DP *)calloc(zd,sizeof(DP));

	utin = (DP *)calloc(zd,sizeof(DP));
	vtin = (DP *)calloc(zd,sizeof(DP));
	wtin = (DP *)calloc(zd,sizeof(DP));
	reys = (DP *)calloc(zd,sizeof(DP));
	dtin = (DP *)calloc(zd,sizeof(DP));

	uv = (DP *)calloc(zd,sizeof(DP));
	uvt = (DP *)calloc(zd,sizeof(DP));
	vw = (DP *)calloc(zd,sizeof(DP));
	vwt = (DP *)calloc(zd,sizeof(DP));

	uavgt = (DP *)calloc(zd,sizeof(DP));
	vavgt = (DP *)calloc(zd,sizeof(DP));
	wavgt = (DP *)calloc(zd,sizeof(DP));
	davgt = (DP *)calloc(zd,sizeof(DP));

	UW = (DP *)calloc(zd,sizeof(DP));
	UWt= (DP *)calloc(zd,sizeof(DP));

	utint = (DP *)calloc(zd,sizeof(DP));
	vtint = (DP *)calloc(zd,sizeof(DP));
	wtint = (DP *)calloc(zd,sizeof(DP));
	reyst = (DP *)calloc(zd,sizeof(DP));
	dtint = (DP *)calloc(zd,sizeof(DP));

	dudx = (DP *)calloc(up,sizeof(DP));
        du = (DP*)calloc(up,sizeof(DP));
        
        sku = (DP *)calloc(zd,sizeof(DP));
        skv = (DP *)calloc(zd,sizeof(DP));
        skw = (DP *)calloc(zd,sizeof(DP));
        skuw = (DP *)calloc(zd,sizeof(DP));
        skut = (DP *)calloc(zd,sizeof(DP));
        skvt = (DP *)calloc(zd,sizeof(DP));
        skwt = (DP *)calloc(zd,sizeof(DP));
        skuwt = (DP *)calloc(zd,sizeof(DP));
        skuvt = (DP *)calloc(zd,sizeof(DP));
        skuv = (DP *)calloc(zd,sizeof(DP));
        skvwt = (DP *)calloc(zd,sizeof(DP));
        skvw = (DP *)calloc(zd,sizeof(DP));
        
        skp = (DP *)calloc(zd,sizeof(DP));
        skpt = (DP *)calloc(zd,sizeof(DP));
        
        kuu = (DP *)calloc(zd,sizeof(DP));
        kuv = (DP *)calloc(zd,sizeof(DP));
        kuw = (DP *)calloc(zd,sizeof(DP));
        kuuw = (DP *)calloc(zd,sizeof(DP));
        kuut = (DP *)calloc(zd,sizeof(DP));
        kuvt = (DP *)calloc(zd,sizeof(DP));
        kuwt = (DP *)calloc(zd,sizeof(DP));
        kuuwt = (DP *)calloc(zd,sizeof(DP));
        kuuv = (DP *)calloc(zd,sizeof(DP));
        kuuvt = (DP *)calloc(zd,sizeof(DP));
        kuvw = (DP *)calloc(zd,sizeof(DP));
        kuvwt = (DP *)calloc(zd,sizeof(DP));
        
        kup = (DP *)calloc(zd,sizeof(DP));
        kupt = (DP *)calloc(zd,sizeof(DP));

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
	printf("In time stats ut=%f nu=%f\n",utauavg,nu);

	sv=fopen("flowrate.txt","r");
        ct=0;
        while (fscanf(sv,"%lf %lf\n",&ttemp,&ubt)!=EOF)
        {
                if ((ttemp>=TSTR)&&(ttemp<=TEND))
                {
                        ubavg+=ubt;
                        ct++;
                }
        }
        fclose(sv);
        ubavg=ubavg/((double)(ct*yd*(zd-2.)));
        printf("ub=%f\n",ubavg);


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
		fread(uavgt,sizeof(DP),zd,sv);
		fread(vavgt,sizeof(DP),zd,sv);
		fread(wavgt,sizeof(DP),zd,sv);
		fread(utint,sizeof(DP),zd,sv);
		fread(vtint,sizeof(DP),zd,sv);
		fread(wtint,sizeof(DP),zd,sv);
		fread(uvt,sizeof(DP),zd,sv);
		fread(reyst,sizeof(DP),zd,sv);
		fread(vwt,sizeof(DP),zd,sv);
		fread(davgt,sizeof(DP),zd,sv);
		fread(dtint,sizeof(DP),zd,sv);
		fread(skut,sizeof(DP),zd,sv);
		fread(skvt,sizeof(DP),zd,sv);
		fread(skwt,sizeof(DP),zd,sv);
		fread(skuvt,sizeof(DP),zd,sv);
		fread(skuwt,sizeof(DP),zd,sv);
		fread(skvwt,sizeof(DP),zd,sv);
		fread(skpt,sizeof(DP),zd,sv);
		fread(kuut,sizeof(DP),zd,sv);
		fread(kuvt,sizeof(DP),zd,sv);
		fread(kuwt,sizeof(DP),zd,sv);
		fread(kuuvt,sizeof(DP),zd,sv);
		fread(kuuwt,sizeof(DP),zd,sv);
		fread(kuvwt,sizeof(DP),zd,sv);
		fread(kupt,sizeof(DP),zd,sv);
		
		
		fread(UWt,sizeof(DP),zd,sv);
		fclose(sv);

		for (a=0;a<zd;a++)
		{
			*(uavg+a) += *(uavgt+a);
			*(vavg+a) += *(vavgt+a);
			*(wavg+a) += *(wavgt+a);
			*(davg+a) += *(davgt+a);

			*(uv+a) += *(uvt+a);
			*(vw+a) += *(vwt+a);

			*(utin+a) += (*(utint+a))*(*(utint+a));
			*(vtin+a) += (*(vtint+a))*(*(vtint+a));
			*(wtin+a) += (*(wtint+a))*(*(wtint+a));
			*(dtin+a) += (*(dtint+a))*(*(dtint+a));
			*(reys+a) += (*(reyst+a));
			*(sku+a) += (*(skut+a));
			*(skv+a) += (*(skvt+a));
			*(skw+a) += (*(skwt+a));
			*(skuv+a) += (*(skuvt+a));
			*(skuw+a) += (*(skuwt+a));
			*(skvw+a) += (*(skvwt+a));
			*(skp+a) += (*(skpt+a));
			*(kuu+a) += (*(kuut+a));
			*(kuv+a) += (*(kuvt+a));
			*(kuw+a) += (*(kuwt+a));
			*(kuuv+a) += (*(kuuvt+a));
			*(kuuw+a) += (*(kuuwt+a));
			*(kuvw+a) += (*(kuvwt+a));
			*(kup+a) += (*(kupt+a));
			*(UW+a) += (*(UWt+a));
		}
	}

	for (k=0;k<zd;k++)
	{
		*(uavg+k) /= num;
		*(vavg+k) /= num;
		*(wavg+k) /= num;
		*(davg+k) /= num;

		*(uv+k) /= num;
		*(vw+k) /= num;
		*(UW+k) /= num;

		*(utin+k) /= num;
		*(vtin+k) /= num;
		*(wtin+k) /= num;
		*(reys+k) /= num;
		*(dtin+k) /= num;
		
		*(sku+k) /= num;
		*(skv+k) /= num;
		*(skw+k) /= num;
		*(skuv+k) /= num;
		*(skuw+k) /= num;
		*(skvw+k) /= num;
		*(skp+k) /= num;
		*(kuu+k) /= num;
		*(kuv+k) /= num;
		*(kuw+k) /= num;
		*(kuuv+k) /= num;
		*(kuuw+k) /= num;
		*(kuvw+k) /= num;
		*(kup+k) /= num;
	}

	for (k=0;k<up;k++)
        *(du+k) = 0.5*((*(uavg+k))+(*(uavg+zd-1-k)));

    *(dudx+0)=(-2.*(*(du+1))+3.*(*(du+2))-(*(du+3)));//((*(z+2)+1.)-(*(z+1)+1.));
    printf("ut=%f tw=%f\n",sqrt(nu*(*(dudx+0))),nu*(*(dudx+0)));
    *(dudx+1)=(-3.*(*(du+1))+4.*(*(du+2))-(*(du+3)))/2.;///(2.*((*(z+2)+1.)-(*(z+1)+1.)));
//        for (k=2;k<(up-1);k++)
//                *(dudx+k)=(*(du+k+1)-*(du+k-1))/2.;//(*(z+k+1)-*(z+k-1));
    for (k=2;k<(up-2);k++)
                *(dudx+k)=(-3.*(*(du+k))+4.*(*(du+k+1))-(*(du+k+2)))/2.;
    *(dudx+up-2)=(*(du+up-2)-*(du+up-3));
    *(dudx+up-1)=(*(du+up-1)-*(du+up-2));//(*(z+up-1)-*(z+up-2));
        
    utitsown=sqrt(nu*(*(dudx+0)));


	sprintf(stri,TEXT);
//	Ret=1.;

	sprintf(fn,"Uavg-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"U<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(uavg+k))+(*(uavg+zd-1-k)))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Uavg-%s-w-its-ut.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"U<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(uavg+k))+(*(uavg+zd-1-k)))/utitsown);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Uavg-Ubulk-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"U/U<sub>bulk</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(uavg+k))/ubavg);
    fprintf(tcp,"\n");
    fclose(tcp);


	sprintf(fn,"Us-b-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
//        fprintf(tcp,"TITLE= \"DNS Results\"\n");
//        fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"U<sup>+</sup>\"\n");
//        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    val=((0.5*((*(uavg+1))+(*(uavg+zd-2))))/utauavg)-(*(z+1)+1.)*Ret;
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(uavg+k))+(*(uavg+zd-1-k)))/utauavg-val);
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"Vavg-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"V<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(vavg+k))+(*(vavg+zd-1-k)))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Vavg-Ubulk-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"V/U<sub>bulk</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(vavg+k))/ubavg);
    fprintf(tcp,"\n");
    fclose(tcp);


	sprintf(fn,"Wavg-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"W<sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(wavg+k))+(*(wavg+zd-1-k)))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Wavg-Ub.dat");
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/h\", \"W<sup>+</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(wavg+k))/ubavg);
    fprintf(tcp,"\n");
    fclose(tcp);

    sprintf(fn,"Wavg-Ub-sym.dat");
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/h\", \"W<sup>+</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*((*(wavg+k))+(*(wavg+zd-1-k)))/ubavg);
    for (k=up-2;k>-1;k--)
        fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),0.5*((*(wavg+k))+(*(wavg+zd-1-k)))/ubavg);
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"Davg-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
//	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>r</greek><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"D AVG\", I=%d, F=POINT\n",up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(davg+k))+(*(davg+zd-1-k)))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"utin-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(utin+k))+(*(utin+zd-1-k))))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"utin-h-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(utin+k))+(*(utin+zd-1-k))))/utauavg);
	for (k=up-2;k>-1;k--)
        fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(utin+k))+(*(utin+zd-1-k))))/utauavg);
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"vtin-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(vtin+k))+(*(vtin+zd-1-k))))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"vtin-h-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(vtin+k))+(*(vtin+zd-1-k))))/utauavg);
	for (k=up-2;k>-1;k--)
        fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(vtin+k))+(*(vtin+zd-1-k))))/utauavg);
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"wtin-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(wtin+k))+(*(wtin+zd-1-k))))/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"wtin-h-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"u<sup>+</sup><sub>rms</sub>,v<sup>+</sup><sub>rms</sub>,w<sup>+</sup><sub>rms</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),sqrt(0.5*((*(wtin+k))+(*(wtin+zd-1-k))))/utauavg);
	for (k=up-2;k>-1;k--)
        fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),sqrt(0.5*((*(wtin+k))+(*(wtin+zd-1-k))))/utauavg);
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"KE-sum-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
//        fprintf(tcp,"TITLE= \"DNS Results\"\n");
//        fprintf(tcp,"VARIABLES = \"Z/H\", \"K<sup>+</sup>\"\n");
//        fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),((*(utin+k))+(*(vtin+k))+(*(wtin+k)))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"uw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<uw><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"uw-no-title-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(reys+k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uw-no-title-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*(*(reys+k)-*(reys+zd-1-k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uw-no-title-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(*(reys+k)-*(reys+zd-1-k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uw-pl-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<uw><sup>+</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"uw-h-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"<uw><sup>+</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),-0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg));
	for (k=up-2;k>-1;k--)
        fprintf(tcp,"%f %.9f\n",2.-((*(z+k))+1),0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"UW-h-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"<UW><sup>+</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(UW+k))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"UW-no-title-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(UW+k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"UW-no-title-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*(*(UW+k)-*(UW+zd-1-k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"UW-no-title-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(*(UW+k)-*(UW+zd-1-k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"UW-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"UW<sup>+</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(-(*(UW+k))+(*(UW+zd-1-k)))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"uv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<uv><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*((*(uv+k))-(*(uv+zd-1-k)))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"uv-no-title-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(uv+k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uv-no-title-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*(*(uv+k)-*(uv+zd-1-k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uv-no-title-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(*(uv+k)-*(uv+zd-1-k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"uv-h-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"<uv><sup>+</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
	    fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(uv+k))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"vw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<vw><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*((*(vw+k))+(*(vw+zd-1-k)))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"vw-no-title-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(vw+k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"vw-no-title-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*(*(vw+k)-*(vw+zd-1-k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"vw-no-title-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"<u<sub>i</sub>u<sub>j</sub>><sup>+</sup>\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(*(vw+k)-*(vw+zd-1-k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"vw-h-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"<vw><sup>+</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(vw+k))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"dtin-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"p<sup>'+</sup><sub>rms</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,sqrt(0.5*((*(dtin+k))+(*(dtin+zd-1-k))))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

    sprintf(fn,"dtin-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"p<sup>'+</sup><sub>rms</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)-1.,sqrt((*(dtin+k)))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);	

	sprintf(fn,"dudx-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,nu*(*(dudx+k))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"dudx-nolog-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),nu*(*(dudx+k))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

    sprintf(fn,"dudx+tw-log-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<u'w'><sup>+</sup>+<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg)+nu*(*(dudx+k))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

    sprintf(fn,"dudx+tw-nolog-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"<u'w'><sup>+</sup>+<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg)+nu*(*(dudx+k))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"dudx+tw-log-new-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"<u'w'><sup>+</sup>+<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up-1);
    for (k=1;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(-(*(UW+k))+(*(UW+zd-1-k)))/(utauavg*utauavg)+0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg)+nu*(*(dudx+k))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"dudx+tw-nolog-new-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z/H\", \"<u'w'><sup>+</sup>+<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up-1);
    for (k=1;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),0.5*(-(*(UW+k))+(*(UW+zd-1-k)))/(utauavg*utauavg)+0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg)+nu*(*(dudx+k))/(utauavg*utauavg));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"skew-u-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(sku+k))+(*(sku+zd-1-k))));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"skew-u-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(sku+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"skew-v-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(skv+k))+(*(skv+zd-1-k))));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"skew-v-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(skv+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"skew-w-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(skw+k))-(*(skw+zd-1-k))));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"skew-w-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(skw+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"skew-p-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(skp+k))+(*(skp+zd-1-k))));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"skew-p-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>u</sub>, S<sub>v</sub>, S<sub>w</sub>, S<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(skp+k)));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"skew-all-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    for (k=0;k<up;k++)
        fprintf(tcp,"%d %f %.9f %.9f %.9f %.9f\n",k,((*(z+k))+1)*Ret,0.5*((*(sku+k))+(*(sku+zd-1-k))),0.5*((*(skw+k))-(*(skw+zd-1-k))),0.5*((*(skv+k))+(*(skv+zd-1-k))),0.5*((*(skp+k))+(*(skp+zd-1-k))));
    fprintf(tcp,"\n");
    fclose(tcp);        

    sprintf(fn,"skew-uw-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>uw</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(skuw+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"skew-uv-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>uv</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(skuv+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"skew-vw-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"S<sub>vw</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(skvw+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-u-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(kuu+k))+(*(kuu+zd-1-k))));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-u-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kuu+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-v-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
            fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(kuv+k))+(*(kuv+zd-1-k))));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-v-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kuv+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-w-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(kuw+k))+(*(kuw+zd-1-k))));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-w-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kuw+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-p-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
    for (k=0;k<up;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(kup+k))+(*(kup+zd-1-k))));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-p-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>u</sub>, F<sub>v</sub>, F<sub>w</sub>, F<sub>p</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kup+k)));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"kurt-all-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    for (k=0;k<up;k++)
        fprintf(tcp,"%d %f %.9f %.9f %.9f %.9f\n",k,((*(z+k))+1)*Ret,0.5*((*(kuu+k))+(*(kuu+zd-1-k))),0.5*((*(kuw+k))+(*(kuw+zd-1-k))),0.5*((*(kuv+k))+(*(kuv+zd-1-k))),0.5*((*(kup+k))+(*(kup+zd-1-k))));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-uw-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>uw</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
        fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kuuw+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-uv-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>uv</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
            fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kuuv+k)));
    fprintf(tcp,"\n");
    fclose(tcp);
        
    sprintf(fn,"kurt-vw-long-%s.dat",SLPTXT);
    tcp=fopen(fn,"w");
    fprintf(tcp,"TITLE= \"DNS Results\"\n");
    fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"F<sub>vw</sub>\"\n");
    fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
    for (k=0;k<zd;k++)
            fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(kuvw+k)));
    fprintf(tcp,"\n");
    fclose(tcp);

	sprintf(fn,"stat-pl-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"Z/H\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\",\"<u'w'><sup>+</sup>+<greek>n</greek>dU/dx/u<sub><greek>t</greek></sub><sup>2</sup>\"\n");
	for (k=0;k<up;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",sqrt(0.5*(*(utin+k)+*(utin+zd-1-k)))/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(*(vtin+k)+*(vtin+zd-1-k)))/utauavg);
		fprintf(tcp,"%f ",sqrt(0.5*(*(wtin+k)+*(wtin+zd-1-k)))/utauavg);
		fprintf(tcp,"%f ",0.5*(*(uv+k)-*(uv+zd-1-k))/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(*(reys+k)-*(reys+zd-1-k))/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(*(vw+k)-*(vw+zd-1-k))/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(*(UW+k)-*(UW+zd-1-k))/(utauavg*utauavg));
		fprintf(tcp,"%f ",nu*(*(dudx+k))/(utauavg*utauavg));
        fprintf(tcp,"%f\n",0.5*(-(*(UW+k))+(*(UW+zd-1-k)))/(utauavg*utauavg)+0.5*(-(*(reys+k))+(*(reys+zd-1-k)))/(utauavg*utauavg)+nu*(*(dudx+k))/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"stat-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"VARIABLES =\"Z+\",\"<u><sub>rms</sub><sup>+</sup>\",\"<v><sub>rms</sub><sup>+</sup>\",\"<w><sub>rms</sub><sup>+</sup>\",\"<uv><sub>rms</sub><sup>+</sup>\"");
	fprintf(tcp,",\"<uw><sub>rms</sub><sup>+</sup>\",\"<vw><sub>rms</sub><sup>+</sup>\",\"<<U><W>><sub>rms</sub><sup>+</sup>\"\n");
	for (k=0;k<zd;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",sqrt(*(utin+k))/utauavg);
		fprintf(tcp,"%f ",sqrt(*(vtin+k))/utauavg);
		fprintf(tcp,"%f ",sqrt(*(wtin+k))/utauavg);
		fprintf(tcp,"%f ",(*(uv+k))/(utauavg*utauavg));
		fprintf(tcp,"%f ",(*(reys+k))/(utauavg*utauavg));
		fprintf(tcp,"%f ",(*(vw+k))/(utauavg*utauavg));
		fprintf(tcp,"%f\n",(*(UW+k))/(utauavg*utauavg));
	}
	fprintf(tcp,"\n");
	fclose(tcp);


	
	free(uavg);
	free(vavg);
	free(wavg);
	free(davg);
	free(utin);
	free(vtin);
	free(wtin);
	free(reys);
	free(dtin);
	free(UW);
	free(UWt);

	free(uv);
	free(uvt);
	free(vw);
	free(vwt);

	free(uavgt);
	free(vavgt);
	free(wavgt);
	free(davgt);
	free(utint);
	free(vtint);
	free(wtint);
	free(reyst);
	free(dtint);
	free(du);
	free(dudx);

	free(z);
	
	free(sku);
	free(skv);
	free(skw);
	free(skuw);
	free(skuv);
	free(skvw);
	free(skut);
	free(skvt);
	free(skwt);
	free(skuwt);
	free(skuvt);
	free(skvwt);
	
	free(skp);
	free(skpt);
	
	free(kuu);
	free(kuv);
	free(kuw);
	free(kuuw);
	free(kuut);
	free(kuvt);
	free(kuwt);
	free(kuuwt);
	free(kuuv);
	free(kuuvt);
	free(kuvwt);
	free(kuvw);
	
	free(kup);
	free(kupt);
}
