#ifndef _PARAMS_H

#define _PARAMS_H

#include "typedefs.h"

struct Params{
	
	int nx = 100;
	int ny = 100;
	int nz = 100;
	int max_steps = 100;
	int output_rate = 10;
	
	Double re = 100.;
	Double tau = 1.;
	Double lid_u = 0.01;
	Double lid_w = 0.0;
	Double lid_mag = 0.0;
	Double nu = 0.0;
	Double tol = 1e-2;
};

#endif
