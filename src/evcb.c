#include <R.h>
#include <math.h>



void empcop_c(double *data, double *u, int *n, int *d, int *lu, double (*result)[*lu] ) 
{
	int in, id, iu;
	int p, h;
	double s;
	for(iu=0; iu < (*lu); iu++) {
		s=0;
		for(in=0; in < *n; in++) {
			p = 1;
			for(id=0; id < *d; id++) {
				h = ( data[in + id * (*n)] ) <= ( u[iu + id * (*lu)] ) ;
				p *= h;
			}
			s += p;
		}
		(*result)[iu] = s / (*n );
	}
}
