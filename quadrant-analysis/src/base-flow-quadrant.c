/************************************************************************
 * Quadrant analysis in base turbulent channel flow                     *
 * Amirreza Rastegari                                                   *
 ************************************************************************/
#include "definitions.h"
#include "fun_defs.h"

int main(int argc, char *argv[])
{
	
	int xd,yd,zd,TM;
	int i,j,k,q,a,b,c,nut;
	int cntr;
	long ss,up;
	POINTER V;
	PDATA pd;
	DP Fx=0.,Fy=0.,Fz=0.;
	DP Gx,Gy,Gz;
	FILE *mysave,*tm;
	FILE *sv;
	char fn[20],fn2[20];
	char *ch;
	MPI_Status status;
	int tst,ttt;
	DP U,Q,ux,dens,mdt;
	DP wup,wdo,dummy,davg,mu;
	DP pgrad,dl,dr,tau;
	DP Reb, dplus,Lplus,ubulk;
	int rank,nproc;

	const double ut=0.0400257;
	const double Ret=216.13878;
	Reb=Rebs;
	ubulk=ub;
	
	zd=zdv+2;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	pd.myrank = rank;
	pd.numproc = nproc;
	
	if (pd.myrank != 0)
		pd.left = pd.myrank - 1;
	else
		pd.left = pd.numproc - 1;
	
	if (pd.myrank != (pd.numproc-1))
		pd.right = pd.myrank + 1;
	else
		pd.right = 0;

	yd=ydv;
	xd=1+(xdv/pd.numproc);

	Fx=CFL*CFL*ut*ut*2./((double)(zd-2));
	tau=(3.*CFL*ut*((double)(zd-2))/(2.*Ret))+0.5;
	tau=1./tau;
	Fy=0.;
	Fz=0.;
	Gx=Fx;
	Gy=0.;
	Gz=0.;

	V.s=(int *)calloc(xd*yd*zd,sizeof(int));
	V.vel=(DP *)calloc((xd+3)*yd*zd*3,sizeof(DP));
	V.den=(DP *)calloc((xd+3)*yd*zd,sizeof(DP));
	if (V.vel==NULL)
	{
		fprintf(stderr,"\nUnable to allocate memory for V.vor: %s\n",strerror(errno));
	}

	i=xd/2;ttt=0;

	V=solid_init(xd,yd,zd,V);
	for (tst=TSTR;tst<=TEND;tst+=DT)
	{
		V=reading(xd, yd, zd, V, pd, tst);
		V=qUADrantA(xd,yd,zd,V,pd,tst);
	}

	free(V.vel);
	free(V.s);
	free(V.den);

	MPI_Finalize();

	return(0);
}


