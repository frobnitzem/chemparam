// J. Am. Chem. Soc., Vol. 114, No. 20, 1992.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// power, B, beta, darad, daeps, lj_cut;
#define NB_POWER (0.25)
#define NB_B     (0.2)
#define NB_BETA  (12.0)
#define NB_DARAD (0.8)
#define NB_DAEPS (0.5)

#define SQR(x) ((x)*(x))
#define POW6(x) (SQR(x)*SQR(x)*SQR(x))

// Note: 1-4 QQ scale = 0.75

// Uses MMFF combination rules to compute VDW energy expression constants
// to produce kJ/mol. energy units
// ai == ai (polarizability)
// Ni == Ni
// ri == Ai
// Gi == Gi
// da == i->da + j->da (0, 2, or 3)
//
// Returns: eps
double compute_lr(double ai, double Ni, double ri, double Gi,
                  double aj, double Nj, double rj, double Gj, int da,
                  double *r_ij) {
	double x, y, z;
	double pow_sq, pow_sqx;
	
	ri *= pow(ai, NB_POWER);
	rj *= pow(aj, NB_POWER);
	
	switch(da) {
	    case 0:		// donor - acceptor pair
		x = ri + rj;	// arithmetic mean
		x *= 0.5;
		
		*r_ij = x * NB_DARAD; // scaled radius
		break;
	    case 2:	// donor - donor pair
	    case 3:	// donor - neither pair
		x = ri + rj;	// arithmetic mean
		x *= 0.5;
		
		*r_ij = x;
		break;
	    default:
		x = ri - rj;		// compute (-1)gamma**2 -> x
		y = ri + rj;
		x /= y;
		x = -SQR(x);
		
		x = exp(NB_BETA*x);	// compute R_IJ
		x = 1.0 - x;
		x *= NB_B;
		x += 1.0;
		x *= y*0.5;
		
		*r_ij = x;
		break;
	    }
	
	z = ai/Ni;	// compute lower part	-> y
	y = aj/Nj;
	y = sqrt(z) + sqrt(y);
	y *= POW6(x);
	
		// Energies in kJ/mol.
	x = 4.184*181.16*Gi*Gj*ai*aj;	// compute upper part ->x
	x /= y;				// eps
	
	if(da == 0) // donor-acceptor pair
            return x * NB_DAEPS;
        return x;
}
