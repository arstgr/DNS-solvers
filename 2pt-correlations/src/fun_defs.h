POINTER slip_ridges_init(int xd, int yd, int zd, POINTER V, PDATA pd);
POINTER reading(int xd, int yd, int zd, POINTER V, PDATA pd, long ts, DP Fx, DP ubulk);
POINTER mean_subtract(int xd,int yd,int zd,POINTER V,PDATA pd,long tst,DP Fx);
POINTER ycorrel(int xd,int yd,int zd,POINTER V,PDATA pd,long tst,DP Fx);
POINTER xcorrel(int xd,int yd,int zd,POINTER V,PDATA pd,long tst,DP Fx);
POINTER xcorrel_slp_nslp(int xd,int yd,int zd,POINTER V,PDATA pd,long tst,DP Fx);
