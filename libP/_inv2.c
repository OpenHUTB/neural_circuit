#define EPS	1.0e-30
double fabs(double);

int _inv2(float (*a)[2], float (*ainv)[2])
                            
/*
 * Invert the 2 by 2 matrix 'a' and leave it in 'ainv'.
 * If it is singular, a non-zero value is returned.
 */
  {	register float *aip;
	float det;

	if ((fabs(det = a[0][0]*a[1][1] - a[0][1]*a[1][0])) < EPS)
		return(1);
	aip = (float *)(ainv);

	*aip++ = a[1][1]/det;
	*aip++ = -a[0][1]/det;
	*aip++ = -a[1][0]/det;
	*aip = a[0][0]/det;
	return(0);
  }
