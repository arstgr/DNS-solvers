#include "definitions.h"

void quadrantA(int xd, int yd, int zd, int st, int en)
{
	int i,j,k,a,c;
	FILE *tcp,*wstr,*sv;
	char fn[100];
	int num=(en-st)+1;
	DP *z,dplus,Ret;
	DP utauavg=0.,utvt,time,nu;
	int cnt=0,ct=0,rd;
	int up;
	double ubulk=ub;
	DP ubavg=0.,ubt=0.,twall,val,tau;
	char stri[60];
	DP uttemp,uttemp1,ttemp;
	int NZ=zd;
	DP quadrant[4][NZ],tquadrant[4][NZ],uw[NZ],uwt[NZ];
	
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
	tau=(3./Rebs)*ubulk*(NZ-2.)+0.5;
	nu=(1./6.)*(2.*tau-1.);
	Ret=utauavg*(NZ-2.)*0.5/nu;
	dplus=Ret/(0.5*(zd-2));
	printf("Quadrant: ut=%f nu=%f dplus=%f Ret=%f\n",utauavg,nu,dplus,Ret);
	
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
	
	for (a=0;a<4;a++)
	{
		for (k=0;k<NZ;k++)
		{
			quadrant[a][k] = 0.;
			tquadrant[a][k] = 0.;
		}
	}
	for (k=0;k<NZ;k++)
	  uw[k] = 0.;
	
	for (a=st;a<(en+1);a+=DT)
	{
		sprintf(fn,"dquadrant.%.3d.%.4d",0,a);	
		sv=fopen(fn,"rb");
		fread(&tquadrant[0][0],sizeof(DP),4*NZ,sv);
		fread(&uwt[0],sizeof(DP),NZ,sv);
		fclose(sv);
		
		for (k=0;k<NZ;k++)
		{
			quadrant[0][k] += tquadrant[0][k];
			quadrant[1][k] += tquadrant[1][k];
			quadrant[2][k] += tquadrant[2][k];
			quadrant[3][k] += tquadrant[3][k];
			
			uw[k] += uwt[k];
		}
	}
	
	for (k=0;k<NZ;k++)
	{
		quadrant[0][k] /= ((double)num);
		quadrant[1][k] /= ((double)num);
		quadrant[2][k] /= ((double)num);
		quadrant[3][k] /= ((double)num);
		
		uw[k] /= ((double)num);
	}

	sprintf(stri,TEXT);
	sprintf(fn,"dquadrant-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"z<sup>+</sup>\", \"z/h\"");
	fprintf(tcp,", \"q<sub>1</sub><sup>+</sup>\", \"q<sub>2</sub><sup>+</sup>\", \"q<sub>3</sub><sup>+</sup>\", \"q<sub>4</sub><sup>+</sup>\"");
	fprintf(tcp,", \"q<sub>1</sub><sup>uw</sup>\", \"q<sub>2</sub><sup>uw</sup>\", \"q<sub>3</sub><sup>uw</sup>\", \"q<sub>4</sub><sup>uw</sup>\"");
	fprintf(tcp,", \"<uw><sup>+</sup>\", \"ub/ut\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up-1);
	fprintf(tcp,"0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. %f\n",ubavg/utauavg);
	for (k=1;k<(up-1);k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));
		
		fprintf(tcp,"%f ",0.5*(quadrant[0][k]-quadrant[3][NZ-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(quadrant[1][k]-quadrant[2][NZ-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(quadrant[2][k]-quadrant[1][NZ-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",0.5*(quadrant[3][k]-quadrant[0][NZ-1-k])/(utauavg*utauavg));
		
		fprintf(tcp,"%f ",0.5*(quadrant[0][k]-quadrant[3][NZ-1-k])/(0.5*(uw[k]-uw[NZ-1-k])));
		fprintf(tcp,"%f ",0.5*(quadrant[1][k]-quadrant[2][NZ-1-k])/(0.5*(uw[k]-uw[NZ-1-k])));
		fprintf(tcp,"%f ",0.5*(quadrant[2][k]-quadrant[1][NZ-1-k])/(0.5*(uw[k]-uw[NZ-1-k])));
		fprintf(tcp,"%f ",0.5*(quadrant[3][k]-quadrant[0][NZ-1-k])/(0.5*(uw[k]-uw[NZ-1-k])));
		
		fprintf(tcp,"%f ",0.5*(uw[k]-uw[NZ-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",ubavg/utauavg);
	}
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"dquadrant-%s-mod.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"z<sup>+</sup>\", \"z/h\"");
	fprintf(tcp,", \"q<sub>1</sub><sup>+</sup>\", \"q<sub>2</sub><sup>+</sup>\", \"q<sub>3</sub><sup>+</sup>\", \"q<sub>4</sub><sup>+</sup>\"");
	fprintf(tcp,", \"q<sub>1</sub><sup>uw</sup>\", \"q<sub>2</sub><sup>uw</sup>\", \"q<sub>3</sub><sup>uw</sup>\", \"q<sub>4</sub><sup>uw</sup>\"");
	fprintf(tcp,", \"<uw><sup>+</sup>\", \"ub/ut\"\n");
	fprintf(tcp,"ZONE T=Base, I=%d, F=POINT\n",up-1);
	fprintf(tcp,"0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. %f\n",ubavg/utauavg);
	for (k=1;k<(up-1);k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(z+k))+1));

		fprintf(tcp,"%f ",(quadrant[0][k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(quadrant[1][k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(quadrant[2][k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(quadrant[3][k])/(utauavg*utauavg));

		fprintf(tcp,"%f ",(quadrant[0][k])/(0.5*(uw[k]-uw[NZ-1-k])));
		fprintf(tcp,"%f ",(quadrant[1][k])/(0.5*(uw[k]-uw[NZ-1-k])));
		fprintf(tcp,"%f ",(quadrant[2][k])/(0.5*(uw[k]-uw[NZ-1-k])));
		fprintf(tcp,"%f ",(quadrant[3][k])/(0.5*(uw[k]-uw[NZ-1-k])));

		fprintf(tcp,"%f ",0.5*(uw[k]-uw[NZ-1-k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",ubavg/utauavg);
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"dquadrant-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	for (k=0;k<NZ;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1));

		fprintf(tcp,"%f ",(quadrant[0][k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(quadrant[1][k])/(utauavg*utauavg));
		fprintf(tcp,"%f ",(quadrant[2][k])/(utauavg*utauavg));
		fprintf(tcp,"%f\n",(quadrant[3][k])/(utauavg*utauavg));
	}
	fprintf(tcp,"%f %f 0. 0. 0. 0. 0. 0. 0. 0. 0. %f\n",((*(z+up-1))+1)*Ret,((*(z+up-1))+1),ubavg/utauavg);
	fprintf(tcp,"\n");
	fclose(tcp);
}
