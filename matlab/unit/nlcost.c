#include "mex.h"

void mexFunction(
		int nlhs,       mxArray *plhs[],
		int nrhs, const mxArray *prhs[]
		) {
	double *y, *u, *theta0, *bb;
	double *cost, e, x;
	int N, t;
 
	y = mxGetPr(prhs[0]);
	N = mxGetM(prhs[0]);
	u = mxGetPr(prhs[1]);
	theta0 = mxGetPr(prhs[2]);
	bb = mxGetPr(prhs[3]);
	
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	cost = mxGetPr(plhs[0]);
	cost[0] = 0.0;
	
	x = 0.0;
	for(t=0;t<N;t++){
		e = y[t]-theta0[3]*x*x;
		cost[0] += e*e;
		x = theta0[0]*x + bb[0]*(x/(1+x*x)) + theta0[2]*u[t];
	}
	
	cost[0] = cost[0] / N;
}