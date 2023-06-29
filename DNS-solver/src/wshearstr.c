#include "definitions.h"

POINTER wshearstr(int xd, int yd, int zd, POINTER V, DP Gx, DP tau, int TM, long cntr, MPI_Comm cart_grid)
{
	int a,mm;
	int k=1,i,j;
	DP stress=0.,ux=0.,dens=0.;
	DP stup=0.,stdown=0.,dummy=0.;
	DP mu,davg=0.;
	FILE *fn;
	char fl[30];
	double dt,Re,umax,dplus;
	double temp,dq;
	double fac,vel,rho,coef=2./3.;
	double ubulk=ub;
	DP ux1,ux2,ux3;
	int in,ip;
	DP dwdx,wp,wn,w,rr,rl,pgrad,pgradt;
	int rank;

	umax=1.;//uts*(2.5*log(Re)+5.5);
	dt=1./((double)(0.5*(zd-2.)/CFL));

	mu=(1./6.)*(2./tau-1);

	for (i=0;i<(xd-1);i++)
	{
		for (j=0;j<(yd-1);j++)
		{
			k=1;
			ux1=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
				ux1+= E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
			}
			ux1 = ux1/dens;
			davg += dens;
			k=2;
			ux2=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
				ux2+= E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
			}
			ux2 = ux2/dens;

			k=3;
			ux3=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
				ux3+= E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
			}
			ux3 = ux3/dens;

			stdown += (-2.*ux1+3.*ux2-ux3);
		}
	}

	davg /= ((double)((xd-1)*(yd-1)*XDIM*YDIM));
	stdown /= ((double)((xd-1)*(yd-1)*XDIM*YDIM));

	MPI_Allreduce(&stdown,&dummy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	stdown=dummy;

	MPI_Allreduce(&davg,&dummy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	davg=dummy;

	stdown *= mu*davg;

	k=zd-2;
	davg=0.;
	for (i=0;i<(xd-1);i++)
	{
		for (j=0;j<(yd-1);j++)
		{
			k=zd-2;
			ux1=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
				ux1+= E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
			}
			ux1 = ux1/dens;
			k=zd-3;
			ux2=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
				ux2+= E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
			}
			ux2 = ux2/dens;

			k=zd-4;
			ux3=0.;dens=0.;
			for (a=0;a<19;a++)
			{
				dens+=Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
				ux3+= E0[a]*Fb(V.ftemp,i,j,k,a,xd,yd,zd,19);
			}
			ux3 = ux3/dens;

			stup += -(2.*ux1-3.*ux2+ux3);
			davg += dens;
		}
	}

	davg /= ((double)((xd-1)*(yd-1)*XDIM*YDIM));
	stup /= ((double)((xd-1)*(yd-1)*XDIM*YDIM));
	dummy=0.;

	MPI_Allreduce(&stup,&dummy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	stup=dummy;

	MPI_Allreduce(&davg,&dummy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	davg=dummy;

	pgrad=V.force;

	pgrad = pgrad*(zd-2)/2.;

	stup *= mu*davg;

	stress=0.5*(stup+stdown);
	if (isnan(stress))
	{
		fprintf(stderr,"Eroor: NAN, found in wall shear stress\n");
		abort();
	}

	MPI_Comm_rank(cart_grid, &rank);
	if (rank==0)
	{
		sprintf(fl,"str-up.txt");
		fn=fopen(fl,"a");
		fprintf(fn,"%f %.12f\n",((double)TM+dt*cntr),stup);
		fclose(fn);

		sprintf(fl,"str-do.txt");
		fn=fopen(fl,"a");
		fprintf(fn,"%f %.12f\n",((double)TM+dt*cntr),stdown);
		fclose(fn);

		sprintf(fl,"w-str.txt");
		fn=fopen(fl,"a");
		fprintf(fn,"%f %.12f\n",((double)TM+dt*cntr),stress);
		fclose(fn);

		sprintf(fl,"p-grad.txt");
		fn=fopen(fl,"a");
		fprintf(fn,"%f %.12f\n",((double)TM+dt*cntr),pgrad);
                fclose(fn);
	}

	return V;
}