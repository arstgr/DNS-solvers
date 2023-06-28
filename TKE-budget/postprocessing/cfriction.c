#include "definitions.h"

void CFriction(int xd, int yd, int zd, int st, int en)
{
	FILE *sv,*tcp;
        int i,j,k,a;
        DP uttemp1,uttemp,utauavg,ttemp;
        int ct=0,ctr,num,mnum;
        char stri[60];
        double ub=0.,cf,cft,ubt,nu,tau,Ret;
        char txt[40];

        sprintf(txt,TEXT);

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
	ctr=0;
	cft=0.;
	num=0;
	sv=fopen("w-str.txt","r");
	while (fscanf(sv,"%lf %lf\n",&ttemp,&uttemp)!=EOF)
        {
                cf=uttemp/(0.5*ub*ub);
		num++;
                if (ttemp>=TSTR)
                {
                        cft+=cf;
                        ctr++;
                }
        }
        cft/=((double)ctr);
	rewind(sv);
	tcp=fopen("Cf-w.dat","w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"T\", \"C<sub>f</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s C<sub>f</sub>-w avg=%f DR=%f\", I=%d, F=POINT\n",txt,cft,(0.007981-cft)*100./0.007981,num);
	while (fscanf(sv,"%lf %lf\n",&ttemp,&uttemp)!=EOF)
        {
                cf=uttemp/(0.5*ub*ub);
                fprintf(tcp,"%lf %.12lf\n",ttemp,cf);
        }
        fclose(sv);
        fclose(tcp);
        printf("average Cf=%f based on wstr\n",cft);
	printf("DR=%f\n",(0.007981-cft)*100./0.007981);


        sv=fopen("p-grad.txt","r");
	ctr=0;
	cft=0.;
	num=0;
	while (fscanf(sv,"%lf %lf\n",&ttemp,&uttemp)!=EOF)
        {
                cf=uttemp/(0.5*ub*ub);
		num++;
                if (ttemp>=TSTR)
                {
                        cft+=cf;
                        ctr++;
                }
        }
        cft/=((double)ctr);
        rewind(sv);
        tcp=fopen("Cf-p.dat","w");
        fprintf(tcp,"TITLE= \"DNS Results\"\n");
        fprintf(tcp,"VARIABLES = \"T\", \"C<sub>f</sub>\"\n");
        fprintf(tcp,"ZONE T=\"%s C<sub>f</sub>-p avg=%f DR=%f\", I=%d, F=POINT\n",txt,cft,(0.007686-cft)*100.0/0.007686,num);
        while (fscanf(sv,"%lf %lf\n",&ttemp,&uttemp)!=EOF)
        {
                cf=uttemp/(0.5*ub*ub);
                fprintf(tcp,"%lf %.12lf\n",ttemp,cf);
        }
        fclose(sv);
        fclose(tcp);
        printf("average Cf=%f based on grad(P)\n",cft);
	printf("DR=%f\n",(0.007686-cft)*100.0/0.007686);

}
