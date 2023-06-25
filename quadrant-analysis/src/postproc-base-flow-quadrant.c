#include "definitions.h"

void quadrantA(int xd, int yd, int zd, int st, int en);

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
	
	quadrantA(xd,yd,zd,TSTR,TEND);
  
  return 0;
}