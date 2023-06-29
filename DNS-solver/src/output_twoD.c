#include "definitions.h"

POINTER output_twoD(int xd, int yd, int zd, POINTER V, DP C, DP dx, DP Fx, DP rhozero, DP tau)
{
	int i=xd/2,j=yd/2,k,a;
	FILE *myfile;
	DP *y;
	DP *jj,dens;
	DP viscosity=C*(2./tau-1.)/6.;

	y=(DP *)calloc(zd,sizeof(DP));
	jj=(DP *)calloc(zd,sizeof(DP));
	
	for (k=0;k<zd;k++)
	  *(y+k) = 0.;
	for (k=1;k<(zd-1);k++)
	{
		*(y+k) = 0.;
		*(y+k)=-6.*ub*((k-0.5)/(zd-2.))*((k-0.5)/(zd-2.)-1);
	}

	for (k=1;k<(zd-1);k++)
	{
		dens=0.;
		for (a=0;a<19;a++)
		{
			dens+=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
			*(jj+k) += E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
		}
	//	*(jj+k) += Fx*0.5;
		*(jj+k) /= dens;
	}
	
	myfile=fopen("2dresults-exact.dat","w");
        fprintf(myfile,"0. %.12f\n",0.);
        for (k=1;k<(zd-1);k++)
                fprintf(myfile,"%.12f %.12lf\n",(k-0.5)/(zd-2.),(*(y+k)));
        fprintf(myfile,"1. 0.\n");
        fclose(myfile);

        myfile=fopen("2dresults-lb.dat","w");
        fprintf(myfile,"0. %.12f\n",*(jj));
        for (k=1;k<(zd-1);k++)
                fprintf(myfile,"%.12f %.12lf\n",(k-0.5)/(zd-2.),(*(jj+k)));
        fprintf(myfile,"1. 0.\n");
        fclose(myfile);
        
	free(y);
	free(jj);
	return V;
}