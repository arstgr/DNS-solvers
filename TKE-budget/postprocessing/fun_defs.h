void Vorticity(int xd, int yd, int zd, int st, int en);
void time_stat(int xd, int yd, int zd, int st, int en);
void Ebudget(int xd, int yd, int zd, int st, int en);
void CFriction(int xd, int yd, int zd, int st, int en);
void uAVG_pl(int xd, int yd, int zd, int st, int en);
void aniso_tr(int xd, int yd, int zd, int st, int en);
void uvfl_pl(int xd, int yd, int zd, int tst);
void aniso_Line_results(int xd, int yd, int zd, int st, int en);
double integral(double *f, int i);
void aniso_slip_noslip_results(int xd, int yd, int zd, int st, int en);
void budget_Line_results(int xd, int yd, int zd, int st, int en);
void budget_Line_results_collected(int xd, int yd, int zd, int st, int en);
void vort_Line_results(int xd, int yd, int zd, int st, int en, int tinst);
void KE_contour(int xd, int yd, int zd, int st, int en);