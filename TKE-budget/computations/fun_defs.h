POINTER solid_init(int xd, int yd, int zd, POINTER V);
POINTER reading(int xd, int yd, int zd, POINTER V, PDATA pd, long ts, DP Fx, DP ubulk);
POINTER uv_fluc(int xd, int yd, int zd, POINTER V, PDATA pd, long ts, DP Fx, int position);
POINTER budget2D(int xd, int yd, int zd, POINTER V, PDATA pd, long ts, DP Fx);
