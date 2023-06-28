#include "definitions.h"

void Ebudget(int xd, int yd, int zd, int st, int en)
{
	int i,j,k,a;
	FILE *sv,*sw;
	FILE *tcp,*wstr;
	char fn[20];
	int num=(en-st)+1;
	DP *z,dplus,Ret;
	DP utauavg=0.,utvt,time,nu;
	int cnt=0,ct=0,rd;
	int up;
	DP ubavg=0.,ubt=0.,twall,val,tau;
	char stri[60];
	DP uttemp,uttemp1,ttemp;
	DP *Pii,*Puu,*Pvv,*Pww,*Puw,*Puv,*Pvw,*Paiuu,*Paivv,*Paiww,*Paiuw;
	DP *Piit,*Puut,*Pvvt,*Pwwt,*Puwt,*Puvt,*Pvwt,*Paiuut,*Paivvt,*Paiwwt,*Paiuwt;
	DP *Paiuv,*Paiuvt,*Paivw,*Paivwt;
	DP *Tiir,*Tuur,*Tvvr,*Twwr,*Tuwr,*Tuvr,*Tvwr;
	DP *Tiirt,*Tuurt,*Tvvrt,*Twwrt,*Tuwrt,*Tuvrt,*Tvwrt;
	DP *Tiipai,*Tuupai,*Tvvpai,*Twwpai,*Tuwpai,*Tuvpai,*Tvwpai;
	DP *Tiipait,*Tuupait,*Tvvpait,*Twwpait,*Tuwpait,*Tuvpait,*Tvwpait;
	DP *Tiis,*Tuus,*Tvvs,*Twws,*Tuws,*Tuvs,*Tvws;
	DP *Tiist,*Tuust,*Tvvst,*Twwst,*Tuwst,*Tuvst,*Tvwst;
	DP *Epii,*Epuu,*Epvv,*Epww,*Epuw,*Epuv,*Epvw;
	DP *Epiit,*Epuut,*Epvvt,*Epwwt,*Epuwt,*Epuvt,*Epvwt;
	DP *KEt,*KE;
	DP *VdKdy,*WdKdz,*Vduudy,*Wduudz,*Vdvvdy,*Wdvvdz,*Vdwwdy,*Wdwwdz,*Vduvdy,*Wduvdz,*Vduwdy,*Wduwdz,*Vdvwdy,*Wdvwdz;
	DP *VdKdyt,*WdKdzt,*Vduudyt,*Wduudzt,*Vdvvdyt,*Wdvvdzt,*Vdwwdyt,*Wdwwdzt,*Vduvdyt,*Wduvdzt,*Vduwdyt,*Wduwdzt,*Vdvwdyt,*Wdvwdzt;
	DP coef;
	DP p1[zd],p2[zd],p3[zd],p4[zd],p5[zd],p6[zd],p7[zd],p8[zd],p9[zd],p10[zd],p11[zd],p12[zd],p13[zd],p14[zd],p15[zd],p16[zd],p17[zd],p18[zd];
	DP p1t[zd],p2t[zd],p3t[zd],p4t[zd],p5t[zd],p6t[zd],p7t[zd],p8t[zd],p9t[zd],p10t[zd],p11t[zd],p12t[zd],p13t[zd],p14t[zd],p15t[zd],p16t[zd],p17t[zd],p18t[zd];
	
	Pii = (DP *)calloc(zd,sizeof(DP));
	Puu = (DP *)calloc(zd,sizeof(DP));
	Pvv = (DP *)calloc(zd,sizeof(DP));
	Pww = (DP *)calloc(zd,sizeof(DP));
	Puw = (DP *)calloc(zd,sizeof(DP));
	Puv = (DP *)calloc(zd,sizeof(DP));
	Pvw = (DP *)calloc(zd,sizeof(DP));
	
	Paiuu = (DP *)calloc(zd,sizeof(DP));
	Paivv = (DP *)calloc(zd,sizeof(DP));
	Paiww = (DP *)calloc(zd,sizeof(DP));
	Paiuw = (DP *)calloc(zd,sizeof(DP));
	Paiuv = (DP *)calloc(zd,sizeof(DP));
	Paivw = (DP *)calloc(zd,sizeof(DP));
	
	Piit = (DP *)calloc(zd,sizeof(DP));
	Puut = (DP *)calloc(zd,sizeof(DP));
	Pvvt = (DP *)calloc(zd,sizeof(DP));
	Pwwt = (DP *)calloc(zd,sizeof(DP));
	Puwt = (DP *)calloc(zd,sizeof(DP));
	Puvt = (DP *)calloc(zd,sizeof(DP));
	Pvwt = (DP *)calloc(zd,sizeof(DP));
	
	Paiuut = (DP *)calloc(zd,sizeof(DP));
	Paivvt = (DP *)calloc(zd,sizeof(DP));
	Paiwwt = (DP *)calloc(zd,sizeof(DP));
	Paiuwt = (DP *)calloc(zd,sizeof(DP));
	Paiuvt = (DP *)calloc(zd,sizeof(DP));
	Paivwt = (DP *)calloc(zd,sizeof(DP));
	
	Tiir = (DP *)calloc(zd,sizeof(DP));
	Tuur = (DP *)calloc(zd,sizeof(DP));
	Tvvr = (DP *)calloc(zd,sizeof(DP));
	Twwr = (DP *)calloc(zd,sizeof(DP));
	Tuwr = (DP *)calloc(zd,sizeof(DP));
	Tuvr = (DP *)calloc(zd,sizeof(DP));
	Tvwr = (DP *)calloc(zd,sizeof(DP));
	
	Tiirt = (DP *)calloc(zd,sizeof(DP));
	Tuurt = (DP *)calloc(zd,sizeof(DP));
	Tvvrt = (DP *)calloc(zd,sizeof(DP));
	Twwrt = (DP *)calloc(zd,sizeof(DP));
	Tuwrt = (DP *)calloc(zd,sizeof(DP));
	Tuvrt = (DP *)calloc(zd,sizeof(DP));
	Tvwrt = (DP *)calloc(zd,sizeof(DP));
	
	Tiipai = (DP *)calloc(zd,sizeof(DP));
	Tuupai = (DP *)calloc(zd,sizeof(DP));
	Tvvpai = (DP *)calloc(zd,sizeof(DP));
	Twwpai = (DP *)calloc(zd,sizeof(DP));
	Tuwpai = (DP *)calloc(zd,sizeof(DP));
	Tuvpai = (DP *)calloc(zd,sizeof(DP));
	Tvwpai = (DP *)calloc(zd,sizeof(DP));
	
	Tiipait = (DP *)calloc(zd,sizeof(DP));
	Tuupait = (DP *)calloc(zd,sizeof(DP));
	Tvvpait = (DP *)calloc(zd,sizeof(DP));
	Twwpait = (DP *)calloc(zd,sizeof(DP));
	Tuwpait = (DP *)calloc(zd,sizeof(DP));
	Tuvpait = (DP *)calloc(zd,sizeof(DP));
	Tvwpait = (DP *)calloc(zd,sizeof(DP));
	
	Tiis = (DP *)calloc(zd,sizeof(DP));
	Tuus = (DP *)calloc(zd,sizeof(DP));
	Tvvs = (DP *)calloc(zd,sizeof(DP));
	Twws = (DP *)calloc(zd,sizeof(DP));
	Tuws = (DP *)calloc(zd,sizeof(DP));
	Tuvs = (DP *)calloc(zd,sizeof(DP));
	Tvws = (DP *)calloc(zd,sizeof(DP));
	
	Tiist = (DP *)calloc(zd,sizeof(DP));
	Tuust = (DP *)calloc(zd,sizeof(DP));
	Tvvst = (DP *)calloc(zd,sizeof(DP));
	Twwst = (DP *)calloc(zd,sizeof(DP));
	Tuwst = (DP *)calloc(zd,sizeof(DP));
	Tuvst = (DP *)calloc(zd,sizeof(DP));
	Tvwst = (DP *)calloc(zd,sizeof(DP));
	
	Epii = (DP *)calloc(zd,sizeof(DP));
	Epuu = (DP *)calloc(zd,sizeof(DP));
	Epvv = (DP *)calloc(zd,sizeof(DP));
	Epww = (DP *)calloc(zd,sizeof(DP));
	Epuw = (DP *)calloc(zd,sizeof(DP));
	Epuv = (DP *)calloc(zd,sizeof(DP));
	Epvw = (DP *)calloc(zd,sizeof(DP));
	
	Epiit = (DP *)calloc(zd,sizeof(DP));
	Epuut = (DP *)calloc(zd,sizeof(DP));
	Epvvt = (DP *)calloc(zd,sizeof(DP));
	Epwwt = (DP *)calloc(zd,sizeof(DP));
	Epuwt = (DP *)calloc(zd,sizeof(DP));
	Epuvt = (DP *)calloc(zd,sizeof(DP));
	Epvwt = (DP *)calloc(zd,sizeof(DP));
	
	KEt = (DP *)calloc(zd,sizeof(DP));
	KE = (DP *)calloc(zd,sizeof(DP));
	
	VdKdy = (DP *)calloc(zd,sizeof(DP));
	WdKdz = (DP *)calloc(zd,sizeof(DP));
	Vduudy = (DP *)calloc(zd,sizeof(DP));
	Wduudz = (DP *)calloc(zd,sizeof(DP));
	Vdvvdy = (DP *)calloc(zd,sizeof(DP));
	Wdvvdz = (DP *)calloc(zd,sizeof(DP));
	Vdwwdy = (DP *)calloc(zd,sizeof(DP));
	Wdwwdz = (DP *)calloc(zd,sizeof(DP));
	Vduvdy = (DP *)calloc(zd,sizeof(DP));
	Wduvdz = (DP *)calloc(zd,sizeof(DP));
	Vduwdy = (DP *)calloc(zd,sizeof(DP));
	Wduwdz = (DP *)calloc(zd,sizeof(DP));
	Vdvwdy = (DP *)calloc(zd,sizeof(DP));
	Wdvwdz = (DP *)calloc(zd,sizeof(DP));
	
	VdKdyt = (DP *)calloc(zd,sizeof(DP));
	WdKdzt = (DP *)calloc(zd,sizeof(DP));
	Vduudyt = (DP *)calloc(zd,sizeof(DP));
	Wduudzt = (DP *)calloc(zd,sizeof(DP));
	Vdvvdyt = (DP *)calloc(zd,sizeof(DP));
	Wdvvdzt = (DP *)calloc(zd,sizeof(DP));
	Vdwwdyt = (DP *)calloc(zd,sizeof(DP));
	Wdwwdzt = (DP *)calloc(zd,sizeof(DP));
	Vduvdyt = (DP *)calloc(zd,sizeof(DP));
	Wduvdzt = (DP *)calloc(zd,sizeof(DP));
	Vduwdyt = (DP *)calloc(zd,sizeof(DP));
	Wduwdzt = (DP *)calloc(zd,sizeof(DP));
	Vdvwdyt = (DP *)calloc(zd,sizeof(DP));
	Wdvwdzt = (DP *)calloc(zd,sizeof(DP));

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
	printf("In Budget ut=%f nu=%f\n",utauavg,nu);

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

	for (a=0;a<zd;a++)
	{
		p1[a] =0.;
		p2[a] =0.;
		p3[a] =0.;
		p4[a] =0.;
		p5[a] =0.;
		p6[a] =0.;
		p7[a] =0.;
		p8[a] =0.;
		p9[a] =0.;
		p10[a] =0.;
		p11[a] =0.;
		p12[a] =0.;
		p13[a] =0.;
		p14[a] =0.;
		p15[a] =0.;
		p16[a] =0.;
		p17[a] =0.;
		p18[a] =0.;
	}
	
	for (k=st;k<(en+1);k+=DT)
	{
		if (SLPF==0)
			sprintf(fn,"KE-Budget.%.3d.%.4d",0,k);
		else if (SLPF==-1)
			sprintf(fn,"KE-Budget-no-slip.%.3d.%.4d",0,k);
		else if (SLPF==-2)
			sprintf(fn,"KE-Budget-slip.%.3d.%.4d",0,k);
		else if (SLPF==-3)
			sprintf(fn,"KE-Budget-2d.%.3d.%.4d",0,k);
			
		sv=fopen(fn,"rb");

		
		fread(Piit,sizeof(DP),zd,sv);
		fread(Puut,sizeof(DP),zd,sv);
		fread(Pvvt,sizeof(DP),zd,sv);
		fread(Pwwt,sizeof(DP),zd,sv);
		fread(Puvt,sizeof(DP),zd,sv);
		fread(Puwt,sizeof(DP),zd,sv);
		fread(Pvwt,sizeof(DP),zd,sv);
		
		fread(Tiirt,sizeof(DP),zd,sv);
		fread(Tuurt,sizeof(DP),zd,sv);
		fread(Tvvrt,sizeof(DP),zd,sv);
		fread(Twwrt,sizeof(DP),zd,sv);
		fread(Tuvrt,sizeof(DP),zd,sv);
		fread(Tuwrt,sizeof(DP),zd,sv);
		fread(Tvwrt,sizeof(DP),zd,sv);
		
		fread(Tiipait,sizeof(DP),zd,sv);
		fread(Tuupait,sizeof(DP),zd,sv);
		fread(Tvvpait,sizeof(DP),zd,sv);
		fread(Twwpait,sizeof(DP),zd,sv);
		fread(Tuvpait,sizeof(DP),zd,sv);
		fread(Tuwpait,sizeof(DP),zd,sv);
		fread(Tvwpait,sizeof(DP),zd,sv);
		
		fread(Paiuut,sizeof(DP),zd,sv);
		fread(Paivvt,sizeof(DP),zd,sv);
		fread(Paiwwt,sizeof(DP),zd,sv);
		fread(Paiuvt,sizeof(DP),zd,sv);
		fread(Paiuwt,sizeof(DP),zd,sv);
		fread(Paivwt,sizeof(DP),zd,sv);
		
		fread(Tiist,sizeof(DP),zd,sv);
		fread(Tuust,sizeof(DP),zd,sv);
		fread(Tvvst,sizeof(DP),zd,sv);
		fread(Twwst,sizeof(DP),zd,sv);
		fread(Tuvst,sizeof(DP),zd,sv);
		fread(Tuwst,sizeof(DP),zd,sv);
		fread(Tvwst,sizeof(DP),zd,sv);
		
		fread(Epiit,sizeof(DP),zd,sv);
		fread(Epuut,sizeof(DP),zd,sv);
		fread(Epvvt,sizeof(DP),zd,sv);
		fread(Epwwt,sizeof(DP),zd,sv);
		fread(Epuvt,sizeof(DP),zd,sv);
		fread(Epuwt,sizeof(DP),zd,sv);
		fread(Epvwt,sizeof(DP),zd,sv);
		
		fread(KEt,sizeof(DP),zd,sv);
		
		fread(VdKdyt,sizeof(DP),zd,sv);
		fread(WdKdzt,sizeof(DP),zd,sv);
		fread(Vduudyt,sizeof(DP),zd,sv);
		fread(Wduudzt,sizeof(DP),zd,sv);
		fread(Vdvvdyt,sizeof(DP),zd,sv);
		fread(Wdvvdzt,sizeof(DP),zd,sv);
		fread(Vdwwdyt,sizeof(DP),zd,sv);
		fread(Wdwwdzt,sizeof(DP),zd,sv);
		fread(Vduvdyt,sizeof(DP),zd,sv);
		fread(Wduvdzt,sizeof(DP),zd,sv);
		fread(Vduwdyt,sizeof(DP),zd,sv);
		fread(Wduwdzt,sizeof(DP),zd,sv);
		fread(Vdvwdyt,sizeof(DP),zd,sv);
		fread(Wdvwdzt,sizeof(DP),zd,sv);
		
		fclose(sv);

		if (SLPF==0)
			sprintf(fn,"TP-Breakdown.%.3d.%.4d",0,k);
		else if (SLPF==-1)
			sprintf(fn,"TP-Breakdown-no-slip.%.3d.%.4d",0,k);
		else if (SLPF==-2)
			sprintf(fn,"TP-Breakdown-slip.%.3d.%.4d",0,k);
		else if (SLPF==-3)
			sprintf(fn,"TP-Breakdown.%.3d.%.4d",0,k);
			
		sv=fopen(fn,"rb");

		fread(&p1t[0],sizeof(DP),zd,sv);
		fread(&p2t[0],sizeof(DP),zd,sv);
		fread(&p3t[0],sizeof(DP),zd,sv);
		fread(&p4t[0],sizeof(DP),zd,sv);
		fread(&p5t[0],sizeof(DP),zd,sv);
		fread(&p6t[0],sizeof(DP),zd,sv);
		fread(&p7t[0],sizeof(DP),zd,sv);
		fread(&p8t[0],sizeof(DP),zd,sv);
		fread(&p9t[0],sizeof(DP),zd,sv);
		fread(&p10t[0],sizeof(DP),zd,sv);
		fread(&p11t[0],sizeof(DP),zd,sv);
		fread(&p12t[0],sizeof(DP),zd,sv);
		fread(&p13t[0],sizeof(DP),zd,sv);
		fread(&p14t[0],sizeof(DP),zd,sv);
		fread(&p15t[0],sizeof(DP),zd,sv);
		fread(&p16t[0],sizeof(DP),zd,sv);
		fread(&p17t[0],sizeof(DP),zd,sv);
		fread(&p18t[0],sizeof(DP),zd,sv);

		fclose(sv);

		for (a=0;a<zd;a++)
		{
			*(Pii+a) += (*(Piit+a));
			*(Puu+a) += (*(Puut+a));
			*(Pvv+a) += (*(Pvvt+a));
			*(Pww+a) += (*(Pwwt+a));
			*(Puv+a) += (*(Puvt+a));
			*(Puw+a) += (*(Puwt+a));
			*(Pvw+a) += (*(Pvwt+a));
			
			*(Tiir+a) += (*(Tiirt+a));
			*(Tuur+a) += (*(Tuurt+a));
			*(Tvvr+a) += (*(Tvvrt+a));
			*(Twwr+a) += (*(Twwrt+a));
			*(Tuvr+a) += (*(Tuvrt+a));
			*(Tuwr+a) += (*(Tuwrt+a));
			*(Tvwr+a) += (*(Tvwrt+a));
			
			*(Tiipai+a) += (*(Tiipait+a));
			*(Tuupai+a) += (*(Tuupait+a));
			*(Tvvpai+a) += (*(Tvvpait+a));
			*(Twwpai+a) += (*(Twwpait+a));
			*(Tuvpai+a) += (*(Tuvpait+a));
			*(Tuwpai+a) += (*(Tuwpait+a));
			*(Tvwpai+a) += (*(Tvwpait+a));
			
			*(Tiis+a) += (*(Tiist+a));
			*(Tuus+a) += (*(Tuust+a));
			*(Tvvs+a) += (*(Tvvst+a));
			*(Twws+a) += (*(Twwst+a));
			*(Tuvs+a) += (*(Tuvst+a));
			*(Tuws+a) += (*(Tuwst+a));
			*(Tvws+a) += (*(Tvwst+a));

			*(Epii+a) += (*(Epiit+a));
			*(Epuu+a) += (*(Epuut+a));
			*(Epvv+a) += (*(Epvvt+a));
			*(Epww+a) += (*(Epwwt+a));
			*(Epuv+a) += (*(Epuvt+a));
			*(Epuw+a) += (*(Epuwt+a));
			*(Epvw+a) += (*(Epvwt+a));

			*(Paiuu+a) += (*(Paiuut+a));
			*(Paivv+a) += (*(Paivvt+a));
			*(Paiww+a) += (*(Paiwwt+a));
			*(Paiuv+a) += (*(Paiuvt+a));
			*(Paiuw+a) += (*(Paiuwt+a));
			*(Paivw+a) += (*(Paivwt+a));
			
			*(KE+a) += (*(KEt+a));
			
			*(VdKdy+a) += (*(VdKdyt+a));
			*(WdKdz+a) += (*(WdKdzt+a));
			*(Vduudy+a) += (*(Vduudyt+a));
			*(Wduudz+a) += (*(Wduudzt+a));
			*(Vdvvdy+a) += (*(Vdvvdyt+a));
			*(Wdvvdz+a) += (*(Wdvvdzt+a));
			*(Vdwwdy+a) += (*(Vdwwdyt+a));
			*(Wdwwdz+a) += (*(Wdwwdzt+a));
			*(Vduvdy+a) += (*(Vduvdyt+a));
			*(Wduvdz+a) += (*(Wduvdzt+a));
			*(Vduwdy+a) += (*(Vduwdyt+a));
			*(Wduwdz+a) += (*(Wduwdzt+a));
			*(Vdvwdy+a) += (*(Vdvwdyt+a));
			*(Wdvwdz+a) += (*(Wdvwdzt+a));

			p1[a] += p1t[a];
			p2[a] += p2t[a];
			p3[a] += p3t[a];
			p4[a] += p4t[a];
			p5[a] += p5t[a];
			p6[a] += p6t[a];
			p7[a] += p7t[a];
			p8[a] += p8t[a];
			p9[a] += p9t[a];
			p10[a] += p10t[a];
			p11[a] += p11t[a];
			p12[a] += p12t[a];
			p13[a] += p13t[a];
			p14[a] += p14t[a];
			p15[a] += p15t[a];
			p16[a] += p16t[a];
			p17[a] += p17t[a];
			p18[a] += p18t[a];
	
		}
	}

	for (a=0;a<zd;a++)
	{
		*(Pii+a) /= num;
		*(Puu+a) /= num;
		*(Pvv+a) /= num;
		*(Pww+a) /= num;
		*(Puv+a) /= num;
		*(Puw+a) /= num;
		*(Pvw+a) /= num;
			
		*(Tiir+a) /= num;
		*(Tuur+a) /= num;
		*(Tvvr+a) /= num;
		*(Twwr+a) /= num;
		*(Tuvr+a) /= num;
		*(Tuwr+a) /= num;
		*(Tvwr+a) /= num;
			
		*(Tiipai+a) /= num;
		*(Tuupai+a) /= num;
		*(Tvvpai+a) /= num;
		*(Twwpai+a) /= num;
		*(Tuvpai+a) /= num;
		*(Tuwpai+a) /= num;
		*(Tvwpai+a) /= num;
			
		*(Tiis+a) /= num;
		*(Tuus+a) /= num;
		*(Tvvs+a) /= num;
		*(Twws+a) /= num;
		*(Tuvs+a) /= num;
		*(Tuws+a) /= num;
		*(Tvws+a) /= num;
			
		*(Epii+a) /= num;
		*(Epuu+a) /= num;
		*(Epvv+a) /= num;
		*(Epww+a) /= num;
		*(Epuv+a) /= num;
		*(Epuw+a) /= num;
		*(Epvw+a) /= num;

		*(Paiuu+a) /= num;
		*(Paivv+a) /= num;
		*(Paiww+a) /= num;
		*(Paiuv+a) /= num;
		*(Paiuw+a) /= num;
		*(Paivw+a) /= num;
		
		*(KE+a) /= num;

		*(VdKdy+a) /= num;
		*(WdKdz+a) /= num;
		*(Vduudy+a) /= num;
		*(Wduudz+a) /= num;
		*(Vdvvdy+a) /= num;
		*(Wdvvdz+a) /= num;
		*(Vdwwdy+a) /= num;
		*(Wdwwdz+a) /= num;
		*(Vduvdy+a) /= num;
		*(Wduvdz+a) /= num;
		*(Vduwdy+a) /= num;
		*(Wduwdz+a) /= num;
		*(Vdvwdy+a) /= num;
		*(Wdvwdz+a) /= num;

		p1[a] /= num;
		p2[a] /= num;
		p3[a] /= num;
		p4[a] /= num;
		p5[a] /= num;
		p6[a] /= num;
		p7[a] /= num;
		p8[a] /= num;
		p9[a] /= num;
		p10[a] /= num;
		p11[a] /= num;
		p12[a] /= num;
		p13[a] /= num;
		p14[a] /= num;
		p15[a] /= num;
		p16[a] /= num;
		p17[a] /= num;
		p18[a] /= num;
	}
	
	sprintf(stri,TEXT);
//	Ret=1.;

//	strcat(stri," P<sub>ii</sub>");
	sprintf(fn,"Pii-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pii+k))+(*(Pii+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pii-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Pii+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," P<sub>uu</sub>");
	sprintf(fn,"Puu-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Puu+k))+(*(Puu+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Puu-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Puu+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," P<sub>vv</sub>");
	sprintf(fn,"Pvv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pvv+k))+(*(Pvv+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pvv-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Pvv+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," P<sub>ww</sub>");
	sprintf(fn,"Pww-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pww+k))+(*(Pww+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pww-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Pww+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Puv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Puv+k))-(*(Puv+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Puv-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Puv+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," P<sub>uw</sub>");
	sprintf(fn,"Puw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Puw+k))-(*(Puw+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Puw-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Puw+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Pvw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Pvw+k))-(*(Pvw+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Pvw-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Pvw+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ii</sub><sup>(r)</sup>");
	sprintf(fn,"Tiir-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiir+k))+(*(Tiir+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tiir-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tiir+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uu</sub><sup>(r)</sup>");
	sprintf(fn,"Tuur-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuur+k))+(*(Tuur+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tuur-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tuur+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>vv</sub><sup>(r)</sup>");
	sprintf(fn,"Tvvr-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvvr+k))+(*(Tvvr+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tvvr-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tvvr+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ww</sub><sup>(r)</sup>");
	sprintf(fn,"Twwr-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Twwr+k))+(*(Twwr+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Twwr-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Twwr+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Tuvr-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuvr+k))-(*(Tuvr+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tuvr-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tuvr+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uw</sub><sup>(r)</sup>");
	sprintf(fn,"Tuwr-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuwr+k))-(*(Tuwr+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tuwr-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tuwr+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Tvwr-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvwr+k))-(*(Tvwr+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Tvwr-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tvwr+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ii</sub><sup>(<greek>p</greek>)</sup>");
	sprintf(fn,"Tiipai-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiipai+k))+(*(Tiipai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tiipai-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tiipai+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uu</sub><sup>(<greek>p</greek>)</sup>");
	sprintf(fn,"Tuupai-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuupai+k))+(*(Tuupai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tuupai-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tuupai+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>vv</sub><sup>(<greek>p</greek>)</sup>");
	sprintf(fn,"Tvvpai-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvvpai+k))+(*(Tvvpai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tvvpai-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tvvpai+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ww</sub><sup>(<greek>p</greek>)</sup>");
	sprintf(fn,"Twwpai-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Twwpai+k))+(*(Twwpai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Twwpai-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Twwpai+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Tuvpai-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuvpai+k))-(*(Tuvpai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tuvpai-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tuvpai+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uw</sub><sup>(<greek>p</greek>)</sup>");
	sprintf(fn,"Tuwpai-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuwpai+k))-(*(Tuwpai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tuwpai-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tuwpai+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Tvwpai-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvwpai+k))-(*(Tvwpai+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tvwpai-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tvwpai+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ii</sub><sup>(s)</sup>");
	sprintf(fn,"Tiis-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiis+k))+(*(Tiis+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tiis-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tiis+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uu</sub><sup>(s)</sup>");
	sprintf(fn,"Tuus-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuus+k))+(*(Tuus+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tuus-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tuus+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>vv</sub><sup>(s)</sup>");
	sprintf(fn,"Tvvs-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvvs+k))+(*(Tvvs+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tvvs-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tvvs+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>ww</sub><sup>(s)</sup>");
	sprintf(fn,"Twws-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Twws+k))+(*(Twws+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Twws-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Twws+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Tuvs-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuvs+k))-(*(Tuvs+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tuvs-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tuvs+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," T<sub>uw</sub><sup>(s)</sup>");
	sprintf(fn,"Tuws-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuws+k))-(*(Tuws+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tuws-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tuws+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Tvws-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvws+k))-(*(Tvws+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Tvws-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Tvws+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>e</greek><sub>ii</sub>");
	sprintf(fn,"Epii-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Epii+k))+(*(Epii+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epii-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Epii+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>e</greek><sub>uu</sub>");
	sprintf(fn,"Epuu-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Epuu+k))+(*(Epuu+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epuu-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Epuu+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>e</greek><sub>vv</sub>");
	sprintf(fn,"Epvv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Epvv+k))+(*(Epvv+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epvv-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Epvv+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>e</greek><sub>ww</sub>");
	sprintf(fn,"Epww-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Epww+k))+(*(Epww+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epww-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Epww+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Epuv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Epuv+k))-(*(Epuv+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epuv-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Epuv+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>e</greek><sub>uw</sub>");
	sprintf(fn,"Epuw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Epuw+k))-(*(Epuw+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epuw-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Epuw+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Epvw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Epvw+k))-(*(Epvw+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Epvw-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Epvw+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>P</greek><sub>uu</sub>");
	sprintf(fn,"Paiuu-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Paiuu+k))+(*(Paiuu+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Paiuu-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Paiuu+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	coef=nu/(utauavg*utauavg*utauavg*utauavg);
	
//	strcat(stri," <greek>P</greek><sub>vv</sub>");
	sprintf(fn,"Paivv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Paivv+k))+(*(Paivv+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Paivv-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Paivv+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>P</greek><sub>ww</sub>");
	sprintf(fn,"Paiww-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Paiww+k))+(*(Paiww+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Paiww-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Paiww+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Paiuv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Paiuv+k))-(*(Paiuv+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Paiuv-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Paiuv+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
//	strcat(stri," <greek>P</greek><sub>uw</sub>");
	sprintf(fn,"Paiuw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Paiuw+k))-(*(Paiuw+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Paiuw-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Paiuw+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Paivw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Paivw+k))-(*(Paivw+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Paivw-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Paivw+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
		
	sprintf(fn,"VdKdy-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdK/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(VdKdy+k))+(*(VdKdy+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"VdKdy-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdK/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(VdKdy+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"WdKdz-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdK/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(WdKdz+k))+(*(WdKdz+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"WdKdz-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdK/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(WdKdz+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Vduudy-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdUU/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Vduudy+k))+(*(Vduudy+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Vduudy-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdUU/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Vduudy+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Wduudz-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdUU/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Wduudz+k))+(*(Wduudz+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Wduudz-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdUU/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Wduudz+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Vdvvdy-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdVV/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Vdvvdy+k))+(*(Vdvvdy+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Vdvvdy-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdVV/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Vdvvdy+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Wdvvdz-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdVV/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Wdvvdz+k))+(*(Wdvvdz+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Wdvvdz-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdVV/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Wdvvdz+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Vdwwdy-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdWW/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Vdwwdy+k))+(*(Vdwwdy+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Vdwwdy-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdWW/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Vdwwdy+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Wdwwdz-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdWW/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Wdwwdz+k))+(*(Wdwwdz+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Wdwwdz-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdWW/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Wdwwdz+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Vduvdy-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdUV/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Vduvdy+k))+(*(Vduvdy+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Vduvdy-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdUV/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Vduvdy+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Wduvdz-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdUV/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Wduvdz+k))+(*(Wduvdz+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Wduvdz-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdUV/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Wduvdz+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Vduwdy-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdUW/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Vduwdy+k))-(*(Vduwdy+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Vduwdy-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdUW/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Vduwdy+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Wduwdz-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdUW/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Wduwdz+k))-(*(Wduwdz+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Wduwdz-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdUW/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Wduwdz+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Vdvwdy-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdVW/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Vdvwdy+k))+(*(Vdvwdy+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Vdvwdy-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s VdVW/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Vdvwdy+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"Wdvwdz-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdVW/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Wdvwdz+k))+(*(Wdvwdz+zd-1-k)))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);	
	
	sprintf(fn,"Wdvwdz-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s WdVW/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(Wdvwdz+k))*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"conv-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"\"Z/H\" \"VdK/dy\" \"WdK/dz\" \"Vduu/dy\" \"Wduu/dz\" \"Vdvv/dy\" \"Wdvv/dz\" \"Vdww/dy\" \"Wdww/dz\" \"Vduv/dy\" \"Wduv/dz\"");
	fprintf(tcp," \"Vduw/dy\" \"Wduw/dz\" \"Vdvw/dy\" \"Wdvw/dz\"\n");
	for (k=0;k<zd;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",(*(VdKdy+k))*coef);
		fprintf(tcp,"%f ",(*(WdKdz+k))*coef);
		fprintf(tcp,"%f ",(*(Vduudy+k))*coef);
		fprintf(tcp,"%f ",(*(Wduudz+k))*coef);
		fprintf(tcp,"%f ",(*(Vdvvdy+k))*coef);
		fprintf(tcp,"%f ",(*(Wdvvdz+k))*coef);
		fprintf(tcp,"%f ",(*(Vdwwdy+k))*coef);
		fprintf(tcp,"%f ",(*(Wdwwdz+k))*coef);
		fprintf(tcp,"%f ",(*(Vduvdy+k))*coef);
		fprintf(tcp,"%f ",(*(Wduvdz+k))*coef);
		fprintf(tcp,"%f ",(*(Vduwdy+k))*coef);
		fprintf(tcp,"%f ",(*(Wduwdz+k))*coef);
		fprintf(tcp,"%f ",(*(Vdvwdy+k))*coef);
		fprintf(tcp,"%f\n",(*(Wdvwdz+k))*coef);
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"conv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"\"Z/H\" \"VdK/dy\" \"WdK/dz\" \"Vduu/dy\" \"Wduu/dz\" \"Vdvv/dy\" \"Wdvv/dz\" \"Vdww/dy\" \"Wdww/dz\" \"Vduv/dy\" \"Wduv/dz\"");
	fprintf(tcp," \"Vduw/dy\" \"Wduw/dz\" \"Vdvw/dy\" \"Wdvw/dz\"\n");
	for (k=0;k<up;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",0.5*((*(VdKdy+k))+(*(VdKdy+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(WdKdz+k))+(*(WdKdz+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Vduudy+k))+(*(Vduudy+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Wduudz+k))+(*(Wduudz+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Vdvvdy+k))+(*(Vdvvdy+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Wdvvdz+k))+(*(Wdvvdz+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Vdwwdy+k))+(*(Vdwwdy+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Wdwwdz+k))+(*(Wdwwdz+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Vduvdy+k))+(*(Vduvdy+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Wduvdz+k))+(*(Wduvdz+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Vduwdy+k))-(*(Vduwdy+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Wduwdz+k))-(*(Wduwdz+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Vdvwdy+k))+(*(Vdvwdy+zd-1-k)))*coef);
		fprintf(tcp,"%f\n",0.5*((*(Wdvwdz+k))+(*(Wdvwdz+zd-1-k)))*coef);
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"KE-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
//	fprintf(tcp,"TITLE= \"DNS Results\"\n");
//	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"K<sup>+</sup>\"\n");
//	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),(*(KE+k))/(utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"testuu-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
//	fprintf(tcp,"TITLE= \"DNS Results\"\n");
//	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
//	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	fprintf(tcp,"Z+ Puu Tuur Tuupai Paiuu Tuus Euu\n");
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",((*(z+k))+1),(*(Puu+k))*coef,(*(Tuur+k))*coef,(*(Tuupai+k))*coef,(*(Paiuu+k))*coef,(*(Tuus+k))*coef,(*(Epuu+k))*coef,((*(Puu+k))+(*(Tuur+k))+(*(Tuupai+k))+(*(Paiuu+k))+(*(Tuus+k))+(*(Epuu+k))-(*(Vduudy+k))-(*(Wduudz+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"testvv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
//	fprintf(tcp,"TITLE= \"DNS Results\"\n");
//	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
//	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	fprintf(tcp,"Z+ Pvv Tvvr Tvvpai Paivv Tvvs Evv\n");
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",((*(z+k))+1),(*(Pvv+k))*coef,(*(Tvvr+k))*coef,(*(Tvvpai+k))*coef,(*(Paivv+k))*coef,(*(Tvvs+k))*coef,(*(Epvv+k))*coef,((*(Pvv+k))+(*(Tvvr+k))+(*(Tvvpai+k))+(*(Paivv+k))+(*(Tvvs+k))+(*(Epvv+k))-(*(Vdvvdy+k))-(*(Wdvvdz+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"testww-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
//	fprintf(tcp,"TITLE= \"DNS Results\"\n");
//	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
//	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	fprintf(tcp,"Z+ Pww Twwr Twwpai Paiww Twws Eww\n");
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",((*(z+k))+1),(*(Pww+k))*coef,(*(Twwr+k))*coef,(*(Twwpai+k))*coef,(*(Paiww+k))*coef,(*(Twws+k))*coef,(*(Epww+k))*coef,((*(Pww+k))+(*(Twwr+k))+(*(Twwpai+k))+(*(Paiww+k))+(*(Twws+k))+(*(Epww+k))-(*(Vdwwdy+k))-(*(Wdwwdz+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);
	
	sprintf(fn,"testii-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
//	fprintf(tcp,"TITLE= \"DNS Results\"\n");
//	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
//	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	fprintf(tcp,"Z+ Pii Tiir Tiipai Tiis Eii\n");
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f %.9f %.9f %.9f %.9f %.9f\n",((*(z+k))+1),(*(Pii+k))*coef,(*(Tiir+k))*coef,(*(Tiipai+k))*coef,(*(Tiis+k))*coef,(*(Epii+k))*coef,((*(Pii+k))+(*(Tiir+k))+(*(Tiipai+k))+(*(Tiis+k))+(*(Epii+k))-(*(VdKdy+k))-(*(WdKdz+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-ii-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),((*(Tiir+k))+(*(Tiis+k))+(*(Tiipai+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-uu-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),((*(Tuur+k))+(*(Tuus+k))+(*(Tuupai+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-vv-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),((*(Tvvr+k))+(*(Tvvs+k))+(*(Tvvpai+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-ww-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),((*(Twwr+k))+(*(Twws+k))+(*(Twwpai+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-uv-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),((*(Tuvr+k))+(*(Tuvs+k))+(*(Tuvpai+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-uw-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),((*(Tuwr+k))+(*(Tuws+k))+(*(Tuwpai+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-vw-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),((*(Tvwr+k))+(*(Tvws+k))+(*(Tvwpai+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-ii-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tiir+k))+(*(Tiir+zd-1-k))+(*(Tiis+k))+(*(Tiis+zd-1-k))+(*(Tiipai+k))+(*(Tiipai+1-k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-uu-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuur+k))+(*(Tuur+zd-1-k))+(*(Tuus+k))+(*(Tuus+zd-1-k))+(*(Tuupai+k))+(*(Tuupai+1-k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-vv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvvr+k))+(*(Tvvr+zd-1-k))+(*(Tvvs+k))+(*(Tvvs+zd-1-k))+(*(Tvvpai+k))+(*(Tvvpai+1-k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-ww-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Twwr+k))+(*(Twwr+zd-1-k))+(*(Twws+k))+(*(Twws+zd-1-k))+(*(Twwpai+k))+(*(Twwpai+1-k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-uv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuvr+k))-(*(Tuvr+zd-1-k))+(*(Tuvs+k))-(*(Tuvs+zd-1-k))+(*(Tuvpai+k))-(*(Tuvpai+1-k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-uw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tuwr+k))-(*(Tuwr+zd-1-k))+(*(Tuws+k))-(*(Tuws+zd-1-k))+(*(Tuwpai+k))-(*(Tuwpai+1-k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Trns-vw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*((*(Tvwr+k))-(*(Tvwr+zd-1-k))+(*(Tvws+k))-(*(Tvws+zd-1-k))+(*(Tvwpai+k))-(*(Tvwpai+1-k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"testuv-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
//	fprintf(tcp,"TITLE= \"DNS Results\"\n");
//	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
//	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	fprintf(tcp,"Z+ Puv Tuvr Tuvpai Paiuv Tuvs Euv\n");
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",((*(z+k))+1),(*(Puv+k))*coef,(*(Tuvr+k))*coef,(*(Tuvpai+k))*coef,(*(Paiuv+k))*coef,(*(Tuvs+k))*coef,(*(Epuv+k))*coef,((*(Puv+k))+(*(Tuvr+k))+(*(Tuvpai+k))+(*(Paiuv+k))+(*(Tuvs+k))+(*(Epuv+k))-(*(Vduvdy+k))-(*(Wduvdz+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"testuw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
//	fprintf(tcp,"TITLE= \"DNS Results\"\n");
//	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
//	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	fprintf(tcp,"Z+ Puw Tuwr Tuwpai Paiuw Tuws Euw\n");
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",((*(z+k))+1),(*(Puw+k))*coef,(*(Tuwr+k))*coef,(*(Tuwpai+k))*coef,(*(Paiuw+k))*coef,(*(Tuws+k))*coef,(*(Epuw+k))*coef,((*(Puw+k))+(*(Tuwr+k))+(*(Tuwpai+k))+(*(Paiuw+k))+(*(Tuws+k))+(*(Epuw+k))-(*(Vduwdy+k))-(*(Wduwdz+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"testvw-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
//	fprintf(tcp,"TITLE= \"DNS Results\"\n");
//	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
//	fprintf(tcp,"ZONE T=\"%s\", I=%d, F=POINT\n",stri,zd);
	fprintf(tcp,"Z+ Pvw Tvwr Tvwpai Paivw Tvws Evw\n");
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n",((*(z+k))+1),(*(Pvw+k))*coef,(*(Tvwr+k))*coef,(*(Tvwpai+k))*coef,(*(Paivw+k))*coef,(*(Tvws+k))*coef,(*(Epvw+k))*coef,((*(Pvw+k))+(*(Tvwr+k))+(*(Tvwpai+k))+(*(Paivw+k))+(*(Tvws+k))+(*(Epvw+k))-(*(Vdvwdy+k))-(*(Wdvwdz+k)))*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Error-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"Z+ ii uu vv ww uv uw vw\n");
	for (k=0;k<zd;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1));
		fprintf(tcp,"%f ",((*(Pii+k))+(*(Tiir+k))+(*(Tiipai+k))+(*(Tiis+k))+(*(Epii+k))-(*(VdKdy+k))-(*(WdKdz+k)))*coef);
		fprintf(tcp,"%f ",((*(Puu+k))+(*(Tuur+k))+(*(Tuupai+k))+(*(Paiuu+k))+(*(Tuus+k))+(*(Epuu+k))-(*(Vduudy+k))-(*(Wduudz+k)))*coef);
		fprintf(tcp,"%f ",((*(Pvv+k))+(*(Tvvr+k))+(*(Tvvpai+k))+(*(Paivv+k))+(*(Tvvs+k))+(*(Epvv+k))-(*(Vdvvdy+k))-(*(Wdvvdz+k)))*coef);
		fprintf(tcp,"%f ",((*(Pww+k))+(*(Twwr+k))+(*(Twwpai+k))+(*(Paiww+k))+(*(Twws+k))+(*(Epww+k))-(*(Vdwwdy+k))-(*(Wdwwdz+k)))*coef);
		fprintf(tcp,"%f ",((*(Puv+k))+(*(Tuvr+k))+(*(Tuvpai+k))+(*(Paiuv+k))+(*(Tuvs+k))+(*(Epuv+k))-(*(Vduvdy+k))-(*(Wduvdz+k)))*coef);
		fprintf(tcp,"%f ",((*(Puw+k))+(*(Tuwr+k))+(*(Tuwpai+k))+(*(Paiuw+k))+(*(Tuws+k))+(*(Epuw+k))-(*(Vduwdy+k))-(*(Wduwdz+k)))*coef);
		fprintf(tcp,"%f\n",((*(Pvw+k))+(*(Tvwr+k))+(*(Tvwpai+k))+(*(Paivw+k))+(*(Tvws+k))+(*(Epvw+k))-(*(Vdvwdy+k))-(*(Wdvwdz+k)))*coef);
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"Error-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"Z+ ii uu vv ww uv uw vw\n");
	for (k=0;k<up;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",((*(Pii+k))+(*(Tiir+k))+(*(Tiipai+k))+(*(Tiis+k))+(*(Epii+k))-(*(VdKdy+k))-(*(WdKdz+k)))*coef);
		fprintf(tcp,"%f ",((*(Puu+k))+(*(Tuur+k))+(*(Tuupai+k))+(*(Paiuu+k))+(*(Tuus+k))+(*(Epuu+k))-(*(Vduudy+k))-(*(Wduudz+k)))*coef);
		fprintf(tcp,"%f ",((*(Pvv+k))+(*(Tvvr+k))+(*(Tvvpai+k))+(*(Paivv+k))+(*(Tvvs+k))+(*(Epvv+k))-(*(Vdvvdy+k))-(*(Wdvvdz+k)))*coef);
		fprintf(tcp,"%f ",((*(Pww+k))+(*(Twwr+k))+(*(Twwpai+k))+(*(Paiww+k))+(*(Twws+k))+(*(Epww+k))-(*(Vdwwdy+k))-(*(Wdwwdz+k)))*coef);
		fprintf(tcp,"%f ",((*(Puv+k))+(*(Tuvr+k))+(*(Tuvpai+k))+(*(Paiuv+k))+(*(Tuvs+k))+(*(Epuv+k))-(*(Vduvdy+k))-(*(Wduvdz+k)))*coef);
		fprintf(tcp,"%f ",((*(Puw+k))+(*(Tuwr+k))+(*(Tuwpai+k))+(*(Paiuw+k))+(*(Tuws+k))+(*(Epuw+k))-(*(Vduwdy+k))-(*(Wduwdz+k)))*coef);
		fprintf(tcp,"%f\n",((*(Pvw+k))+(*(Tvwr+k))+(*(Tvwpai+k))+(*(Paivw+k))+(*(Tvws+k))+(*(Epvw+k))-(*(Vdvwdy+k))-(*(Wdvwdz+k)))*coef);
	}
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uvdU/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p1[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uwdU/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p2[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vvdV/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p3[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vwdV/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p4[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -wvdW/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p5[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -wwdW/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p6[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P7-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uvdW/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p7[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P8-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uwdW/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p8[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P9-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vwdU/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p9[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P10-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -wwdU/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p10[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P11-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uvdV/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p11[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P12-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uwdV/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p12[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P13-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vvdU/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p13[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P14-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vwdU/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p14[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P15-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vvdW/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p15[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P16-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vwdW/dz\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p16[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P17-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vwdV/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p17[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P18-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z/H\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -wwdV/dy\", I=%d, F=POINT\n",stri,zd);
	for (k=0;k<zd;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1),p18[k]*nu/(utauavg*utauavg*utauavg*utauavg));
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P1-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uvdU/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p1[k]+p1[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P2-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uwdU/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p2[k]+p2[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P3-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vvdV/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p3[k]+p3[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P4-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vwdV/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p4[k]+p4[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P5-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -wvdW/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p5[k]+p5[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P6-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -wwdW/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p6[k]+p6[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P7-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uvdW/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p7[k]-p7[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P8-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uwdW/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p8[k]-p8[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P9-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vwdU/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p9[k]-p9[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P10-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -wwdU/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p10[k]-p10[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P11-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uvdV/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p11[k]-p11[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P12-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -uwdV/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p12[k]+p12[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P13-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vvdU/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p13[k]+p13[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P14-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vwdU/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p14[k]-p14[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P15-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vvdW/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p15[k]-p15[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P16-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vwdW/dz\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p16[k]+p16[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P17-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -vwdV/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p17[k]+p17[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P18-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"TITLE= \"DNS Results\"\n");
	fprintf(tcp,"VARIABLES = \"Z<sup>+</sup>\", \"Energy Budget\"\n");
	fprintf(tcp,"ZONE T=\"%s -wwdV/dy\", I=%d, F=POINT\n",stri,up);
	for (k=0;k<up;k++)
		fprintf(tcp,"%f %.9f\n",((*(z+k))+1)*Ret,0.5*(p18[k]+p18[zd-k-1])*coef);
	fprintf(tcp,"\n");
	fclose(tcp);

	sprintf(fn,"P-B-h-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"\"Z/H\" \"-uvdU/dy\" \"-uwdUdz\" \"-vvdV/dy\" \"-vwdV/dz\" \"-wvdW/dy\" \"-wwdW/dz\" \"-uvdW/dy\" \"-uwdW/dz\" \"-vwdU/dy\" \"-wwdU/dz\"");
	fprintf(tcp," \"-uvdVdy\" \"-uwdV/dz\" \"-vvdU/dy\" \"-vwdU/dz\" \"-vvdW/dy\" \"-vwdW/dz\" \"-vwdV/dy\" \"-wwdV/dz\" \"Pii\" \"Puu\" \"Pvv\" \"Pww\"");
	fprintf(tcp," \"Puv\" \"Puw\" \"Pvw\"\n");
	for (k=0;k<zd;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1));//1
		fprintf(tcp,"%f ",p1[k]*coef);//2
		fprintf(tcp,"%f ",p2[k]*coef);//3
		fprintf(tcp,"%f ",p3[k]*coef);//4
		fprintf(tcp,"%f ",p4[k]*coef);//5
		fprintf(tcp,"%f ",p5[k]*coef);//6
		fprintf(tcp,"%f ",p6[k]*coef);//7
		fprintf(tcp,"%f ",p7[k]*coef);//8
		fprintf(tcp,"%f ",p8[k]*coef);//9
		fprintf(tcp,"%f ",p9[k]*coef);//10
		fprintf(tcp,"%f ",p10[k]*coef);//11
		fprintf(tcp,"%f ",p11[k]*coef);//12
		fprintf(tcp,"%f ",p12[k]*coef);//13
		fprintf(tcp,"%f ",p13[k]*coef);//14
		fprintf(tcp,"%f ",p14[k]*coef);//15
		fprintf(tcp,"%f ",p15[k]*coef);//16
		fprintf(tcp,"%f ",p16[k]*coef);//17
		fprintf(tcp,"%f ",p17[k]*coef);//18
		fprintf(tcp,"%f ",p18[k]*coef);//19
		fprintf(tcp,"%f ",(*(Pii+k))*coef);//20
		fprintf(tcp,"%f ",(*(Puu+k))*coef);//21
		fprintf(tcp,"%f ",(*(Pvv+k))*coef);//22
		fprintf(tcp,"%f ",(*(Pww+k))*coef);//23
		fprintf(tcp,"%f ",(*(Puv+k))*coef);//24
		fprintf(tcp,"%f ",(*(Puw+k))*coef);//25
		fprintf(tcp,"%f\n",(*(Pvw+k))*coef);//26
	}
	fprintf(tcp,"\n");
	fclose(tcp);
	

	sprintf(fn,"P-B-%s.dat",SLPTXT);
	tcp=fopen(fn,"w");
	fprintf(tcp,"\"Z+\" \"-uvdU/dy\" \"-uwdUdz\" \"-vvdV/dy\" \"-vwdV/dz\" \"-wvdW/dy\" \"-wwdW/dz\" \"-uvdW/dy\" \"-uwdW/dz\" \"-vwdU/dy\" \"-wwdU/dz\"");
	fprintf(tcp," \"-uvdVdy\" \"-uwdV/dz\" \"-vvdU/dy\" \"-vwdU/dz\" \"-vvdW/dy\" \"-vwdW/dz\" \"-vwdV/dy\" \"-wwdV/dz\" \"Pii\" \"Puu\" \"Pvv\" \"Pww\"");
	fprintf(tcp," \"Puv\" \"Puw\" \"Pvw\"\n");
	for (k=0;k<up;k++)
	{
		fprintf(tcp,"%f ",((*(z+k))+1)*Ret);
		fprintf(tcp,"%f ",0.5*(p1[k]+p1[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p2[k]+p2[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p3[k]+p3[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p4[k]+p4[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p5[k]+p5[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p6[k]+p6[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p7[k]-p7[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p8[k]-p8[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p9[k]-p9[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p10[k]-p10[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p11[k]-p11[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p12[k]+p12[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p13[k]+p13[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p14[k]-p14[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p15[k]-p15[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p16[k]+p16[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p17[k]+p17[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*(p18[k]+p18[zd-k-1])*coef);
		fprintf(tcp,"%f ",0.5*((*(Pii+k))+(*(Pii+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Puu+k))+(*(Puu+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Pvv+k))+(*(Pvv+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Pww+k))+(*(Pww+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Puv+k))-(*(Puv+zd-1-k)))*coef);
		fprintf(tcp,"%f ",0.5*((*(Puw+k))-(*(Puw+zd-1-k)))*coef);
		fprintf(tcp,"%f\n",0.5*((*(Pvw+k))-(*(Pvw+zd-1-k)))*coef);
	}
	fprintf(tcp,"\n");
	fclose(tcp);
	
	free(Pii);
	free(Puu);
	free(Pvv);
	free(Pww);
	free(Puw);
	free(Puv);
	free(Pvw);
	
	free(Paiuu);
	free(Paivv);
	free(Paiww);
	free(Paiuw);
	free(Paiuv);
	free(Paivw);
	
	free(Piit);
	free(Puut);
	free(Pvvt);
	free(Pwwt);
	free(Puwt);
	free(Puvt);
	free(Pvwt);
	
	free(Paiuut);
	free(Paivvt);
	free(Paiwwt);
	free(Paiuwt);
	free(Paiuvt);
	free(Paivwt);
	
	free(Tiir);
	free(Tuur);
	free(Tvvr);
	free(Twwr);
	free(Tuwr);
	free(Tuvr);
	free(Tvwr);
	
	free(Tiirt);
	free(Tuurt);
	free(Tvvrt);
	free(Twwrt);
	free(Tuwrt);
	free(Tuvrt);
	free(Tvwrt);
	
	free(Tiipai);
	free(Tuupai);
	free(Tvvpai);
	free(Twwpai);
	free(Tuwpai);
	free(Tuvpai);
	free(Tvwpai);
	
	free(Tiipait);
	free(Tuupait);
	free(Tvvpait);
	free(Twwpait);
	free(Tuwpait);
	free(Tuvpait);
	free(Tvwpait);
	
	free(Tiis);
	free(Tuus);
	free(Tvvs);
	free(Twws);
	free(Tuws);
	free(Tuvs);
	free(Tvws);	
	
	free(Tiist);
	free(Tuust);
	free(Tvvst);
	free(Twwst);
	free(Tuwst);
	free(Tuvst);
	free(Tvwst);
	
	free(Epii);
	free(Epuu);
	free(Epvv);
	free(Epww);
	free(Epuw);
	free(Epuv);
	free(Epvw);
	
	free(Epiit);
	free(Epuut);
	free(Epvvt);
	free(Epwwt);
	free(Epuwt);
	free(Epuvt);
	free(Epvwt);
	
	free(VdKdy);
	free(WdKdz);
	free(Vduudy);
	free(Wduudz);
	free(Vdvvdy);
	free(Wdvvdz);
	free(Vdwwdy);
	free(Wdwwdz);
	free(Vduvdy);
	free(Wduvdz);
	free(Vduwdy);
	free(Wduwdz);
	free(Vdvwdy);
	free(Wdvwdz);
	
	free(VdKdyt);
	free(WdKdzt);
	free(Vduudyt);
	free(Wduudzt);
	free(Vdvvdyt);
	free(Wdvvdzt);
	free(Vdwwdyt);
	free(Wdwwdzt);
	free(Vduvdyt);
	free(Wduvdzt);
	free(Vduwdyt);
	free(Wduwdzt);
	free(Vdvwdyt);
	free(Wdvwdzt);
	
	free(z);
}