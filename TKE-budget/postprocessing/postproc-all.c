/************************************************************************
* Amirreza Rastegari                                                    *
* This should be used along with the budget-2d-revised                  *
* This code gives all of the statistics and budgets                     *
*************************************************************************/
#include "definitions.h"
#include "fun_defs.h"

int main()
{
	int xd,yd,zd,TM;
	int i,j,k,q,a,b,c,nut;
	int tagl=1,tagr=2;
	int cntr;
	int ss,up;
	FILE *epp;
	const DP C=1.;
	DP dx=1.;
	DP Fx=0.,Fy=0.,Fz=0.;
	DP Gx,Gy,Gz;
	DP tau=0.;
	DP rhozero=1.,rhozeroinv;
	DP t1,t2,t3,t4,t5;
	char ch;
	FILE *mysave,*tm;
	FILE *sv;
	char fn[20],fn2[20];
	DP utime,umax,dt;
	DP utavg;

	zd=zdv+2;
	xd=xdv;
	yd=ydv;

/* turbulence statistics */
	time_stat(xd,yd,zd,TSTR,TEND);
/* time history of the skin friction coefficient */
//	CFriction(xd,yd,zd,TSTR,TEND);
/* rms vorticity fluctuations */
//	Vorticity(xd,yd,zd,TSTR,TEND);
/* contours of the mean UV velocities at given z planes */
//	uAVG_pl(xd,yd,zd,TSTR,TEND);
/* TKE budget, all terms */
//	Ebudget(xd,yd,zd,TSTR,TEND);
/* contours of fluctuating uv velocities at given z planes */
//	uvfl_pl(xd, yd, zd, 800);
/* anisotropy invariants of the Reynolds shear stress tensor */
//	aniso_tr(xd,yd,zd,TSTR,TEND);
/* anisotropy invariants of the Reynolds shear stress tensor at given y planes */
//	aniso_Line_results(xd,yd,zd,TSTR,TEND);
/* TKE budget at given y planes */
//	budget_Line_results(xd,yd,zd,TSTR,TEND);
/* TKE budget averaged over the periodic patterns of surface micro-texture */
//	budget_Line_results_collected(xd,yd,zd,TSTR,TEND);
/* anisotropy invariants of the Reynolds shear stresses averaged over the periodic patterns of surface micro-texture */
//	aniso_slip_noslip_results(xd,yd,zd,TSTR,TEND);
/* vorticity fluctuations averaged over the periodic patterns of surface micro-texture */
//	vort_Line_results(xd,yd,zd,TSTR,TEND,TINST);
/* contour plots of KE */
//	KE_contour(xd,yd,zd,TSTR,TEND);

	return(0);
}
