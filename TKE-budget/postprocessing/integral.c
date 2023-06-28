#include "definitions.h"

/* Integrates f between 0 to i */
double integral(double *f, int i)
{
	double integral=0.;
	int j,k;

	for (j=2;j<(i-1);j++)
		integral += 2.*f[j];

	integral += f[1];
	integral += f[i];

	integral *= 0.5;

	if (i==1)
		integral = 0.;
//	integral += 0.5*0.5*(f[0]+f[1]);

	return integral;
}
