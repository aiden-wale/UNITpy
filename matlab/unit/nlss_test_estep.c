//mex -O nlss_test_estep.c

#define varint int

#include "mex.h"
#include "math.h"
#include <stdlib.h>

typedef struct {
    varint N, M, n, ny, nu, t;
    double *y,*u,*theta;
} parameters;


/*-------------------------------------------------------------------------*/
//Generate a uniform random number between 0 and 1
double ranf(){
    
    return ((double) rand())/((double)RAND_MAX);
}
/*-------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------*/
//Generate a Gaussian random number with unit variance
double randn(){
    
    double x1, x2, w;
 
    do {
        x1 = 2.0 * ranf() - 1.0;
        x2 = 2.0 * ranf() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    return x1 * w;
}
/*-------------------------------------------------------------------------*/



/*-------------------------------------------------------------------------*/
/*                         USER DEFINED FUNCTIONS START                    */
/*-------------------------------------------------------------------------*/

void initialstates(parameters *pars, double *x){
    varint i;
    
    for(i=0;i<pars->M;i++){
        x[i] = randn();
    }
}


void likelihood(parameters *pars, double *x, double *q){
    varint i;
    double tmp, *th, y;
    th = pars->theta;
    y  = pars->y[pars->t];
    
    for(i=0;i<pars->M;i++){
        tmp  = (y - th[3]*x[i]*x[i])/th[5];
        q[i] = exp(-0.5*tmp*tmp);
    }
}

void prediction(parameters *pars, double *xt, double *x){
    varint i;
    double tmp, *th, u;
    th = pars->theta;
    u  = pars->u[pars->t];
    
    for(i=0;i<pars->M;i++){
        x[i]  = th[0]*xt[i] + th[1]*(xt[i]/(1.0+xt[i]*xt[i])) + th[2]*u + th[4]*randn();
    }
}


double transitionratio(parameters *pars, double *xt, double *xt1){
    varint i;
    double tmp, *th, u;
    th = pars->theta;
    u  = pars->u[pars->t];
    
    tmp = (xt1[0] - th[0]*xt[0] - th[1]*(xt[0]/(1.0+xt[0]*xt[0])) - th[2]*u)/th[4];
    return exp(-0.5*tmp*tmp);
}

/*-------------------------------------------------------------------------*/
/*                         USER DEFINED FUNCTIONS END                      */
/*-------------------------------------------------------------------------*/




/*-------------------------------------------------------------------------*/
//Add a matrix to a structure
double * add2d(mxArray *structure, char *name, varint size1, varint size2) 
{
	double *ptr;
	mxArray *tmp;
	
	
	tmp = mxGetField(structure, 0, name);
	if (tmp!=NULL){
		ptr = mxGetPr(tmp);
		if((size1!=mxGetM(tmp)) || (size2!=mxGetN(tmp))){
			mxDestroyArray(tmp);
			tmp = mxCreateDoubleMatrix(size1, size2, mxREAL);
			ptr = mxGetPr(tmp); mxSetField(structure, 0, name, tmp);
		}
	}
	else {
		mxAddField(structure, name);
		tmp = mxCreateDoubleMatrix(size1, size2, mxREAL);
		ptr = mxGetPr(tmp); mxSetField(structure, 0, name, tmp);
	}
	
	return ptr;
	
}
/*-------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------*/
//Add a 3D matrix to a structure
double * add3d(mxArray *structure, char *name, varint size1, varint size2, varint size3) 
{
	double *ptr;
	mwSize dims[3], numd, *dimd;    
	mxArray *tmp;
	
	tmp = mxGetField(structure, 0, name);
	if (tmp!=NULL){
		//mexPrintf("Object exists\n");
		ptr = mxGetPr(tmp); 
		numd  = mxGetNumberOfDimensions(tmp);
		if (numd==3) {
			dimd = mxGetDimensions(tmp);
			if((size1!=dimd[0]) || (size2!=dimd[1]) || (size3!=dimd[2])){
				mxDestroyArray(tmp);
				dims[0] = size1; dims[1] = size2; dims[2] = size3;
				tmp = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
				ptr  = mxGetPr(tmp); mxSetField(structure, 0, name, tmp);
				//mexPrintf("Object destroyed 1\n");
			}
		} else {
			mxDestroyArray(tmp);
			dims[0] = size1; dims[1] = size2; dims[2] = size3;
			tmp = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
			ptr  = mxGetPr(tmp); mxSetField(structure, 0, name, tmp);
			//mexPrintf("Object destroyed 2\n");
		}
	}
	else {
		mxAddField(structure, name);
		dims[0] = size1; dims[1] = size2; dims[2] = size3;
		tmp = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
		ptr  = mxGetPr(tmp); mxSetField(structure, 0, name, tmp);
		//mexPrintf("Object created\n");
	}
	
	return ptr;
	
}
/*-------------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	varint i, j, k, l, t, n, nu, ny, N, M, tmpi, numtries, numsuccessful, kstore, ii, idxcount, whilecount, iszero;
    varint *idx, *idx2;
	double *theta, *x, *xf, *xs, *u, *y, *q, *epsil, *xt, *xt1, *dtmp;
	double sumq, cu, tmp;
	
    mxArray *mxtmp, *nlss, *Z, *G, *O;
    int fcnlhs, fcnrhs; 
    mxArray *fcplhs, *fcprhs[3];
    
    parameters pars; //Structure of useful parameters
    
    if (nrhs!=3){
        mexErrMsgTxt("This routine requires 3 inputs, Z, M, and O");
    }
    
    //First thing is to duplicate M into G
    plhs[0] = mxDuplicateArray(prhs[1]);
    G       = plhs[0];
    Z       = prhs[0];
    O       = prhs[2];
    
    //Get the number of states i.e. G.nx
	mxtmp  = mxGetField(G, 0, "nx");
	if (mxtmp!=NULL){
		n  = (varint)(mxGetPr(mxtmp)[0]);
	}
	else {
		mexErrMsgTxt("No such field: nx, in the model structure");
	}
    
    //Get the number of inputs i.e. Z.nu
	mxtmp  = mxGetField(Z, 0, "nu");
	if (mxtmp!=NULL){
		nu  = (varint)(mxGetPr(mxtmp)[0]);
	}
	else {
		mexErrMsgTxt("No such field: nu, in the data structure");
	}
    
    
    //Get a pointer to the inputs i.e. Z.u
	mxtmp  = mxGetField(Z, 0, "u");
	if (mxtmp!=NULL){
		u  = mxGetPr(mxtmp);
	}
	else {
		mexErrMsgTxt("No such field: u, in the data structure");
	}
    
    //Get the number of outputs i.e. Z.ny
	mxtmp  = mxGetField(Z, 0, "ny");
	if (mxtmp!=NULL){
		ny  = (varint)(mxGetPr(mxtmp)[0]);
	}
	else {
		mexErrMsgTxt("No such field: nu, in the data structure");
	}
    
    //Get a pointer to the outputs i.e. Z.y
	mxtmp  = mxGetField(Z, 0, "y");
	if (mxtmp!=NULL){
		y  = mxGetPr(mxtmp);
	}
	else {
		mexErrMsgTxt("No such field: y, in the data structure");
	}
    
    
    //Get the number of data points Z.Ny
	mxtmp  = mxGetField(Z, 0, "Ny");
	if (mxtmp!=NULL){
		N  = (varint)(mxGetPr(mxtmp)[0]);
	}
	else {
		mexErrMsgTxt("No such field: Ny, in the data structure");
	}
    
    //Get the number of particles OPT.pnum
	mxtmp  = mxGetField(O, 0, "pnum");
	if (mxtmp!=NULL){
		M = (varint)(mxGetPr(mxtmp)[0]);
	}
	else {
		mexErrMsgTxt("No such field: pnum, in the options structure");
	}
    
    //Get the theta vector M.theta
	mxtmp  = mxGetField(G, 0, "theta");
	if (mxtmp!=NULL){
		theta = mxGetPr(mxtmp);
	}
	else {
		mexErrMsgTxt("No such field: theta, in the model structure");
	}
	
    //Add the filtered and smoothed particles
	xf = add3d(G, "xf", n, M, N+1);
    xs = add3d(G, "xs", n, M, N+1);
    
	
	//Make some room
	x      = (double *)malloc(n*M*sizeof(double));
	q      = (double *)malloc(M*sizeof(double));
	epsil  = (double *)malloc(M*sizeof(double));
	idx    = (varint *)malloc(M*sizeof(varint));
	idx2   = (varint *)malloc(M*sizeof(varint));
    
    //Set the common parameters to pass to all functions
    pars.N     = N;
    pars.M     = M;
    pars.n     = n;
    pars.nu    = nu;
    pars.ny    = ny;
    pars.t     = 0;
    pars.y     = y;
    pars.u     = u;
    pars.theta = theta;
    
	//Call function to set initial states
    initialstates(&pars,x);
    
	//Run main loop
	for(t=0;t<N;t++){
        //Update the time variable in pars
        pars.t = t;
        
		//Call the likelihood function
        likelihood(&pars,x,q);
        
        //Add the likelihoods and normalise weights
        sumq = 0.0;
		for(i=0;i<M;i++){
			sumq = sumq + q[i];
		}
		sumq = 1.0/sumq;
		for(i=0;i<M;i++){
            q[i] = q[i]*sumq;
        }
		
		
		//Set pointer to filtered particles at time t
		xt = &xf[n*(t*M)];
		
		//Resample
		for(i=1;i<M;i++){
            q[i]=q[i]+q[i-1];
        }
		sumq=ranf();
		for(i=0;i<M;i++){
            epsil[i]=(i+sumq)/M;
        }
		i=0;
		for(k=0;k<M;k++){
            while(q[i]<epsil[k]){
                i++;
            }
            for(j=0;j<n;j++){
                xt[n*k+j] = x[n*i+j];
            }
        }
        
		//Do the prediction
        prediction(&pars,xt,x);
	}
	
    //Save the final predicted state into the (N+1)'th time slot of filtered and smoothed states
	for(i=0;i<M;i++){
		for(j=0;j<n;j++){
			xf[n*(i+N*M)+j]   = x[n*i+j];
			xs[n*(i+N*M)+j]   = x[n*i+j];
		}
	}
 	
    //Save the N'th filtered state into the N'th smoothed state
	t=N-1;
	for(i=0;i<M;i++){
		for(j=0;j<n;j++){
			xs[n*(i+t*M)+j] = xf[n*(i+t*M)+j];
		}
	}
	
    //Initialise the idx2 to range from 0 to M-1
	for(i=0;i<M;i++){
		idx2[i]=i;
	}
    
    
	//Run smoother
	for(t=N-2;t>=0;t--) {
        //Update the time variable in pars
        pars.t = t;
        
        //Set some pointers to make indexing easier
		xt1    = &xs[n*((t+1)*M)];
		xt     = &xf[n*t*M];
	
        //For each t we reset idx to -1 which indicates that we have not calculated a likelihood for this index
		for(i=0;i<M;i++){
			idx[i]=-1;
		}
		
        //Loop over each SMOOTHED particle xs(:,j,t) at time t+1
		for(j=0;j<M;j++){
			numsuccessful = 0;
			idxcount=M;
			whilecount=0;
			while(numsuccessful < 1){
				whilecount++;
				ii  = (varint)(idxcount*((double) rand())/((double)RAND_MAX));
				k   = idx2[ii];
				
                //If we have already calculated the likelihood then don't repeat it
  				if(idx[k]==j){
  					tmp    = q[k];
  				} else { //Otherwise calculate the likelihood as a fraction of the max possible value that the PDF can have
					tmp = transitionratio(&pars,&xt[n*k],&xt1[n*j]);
					if(tmp < 1e-16) {
                        //If the likelihood is zero then record this as the particle will never be chosen
						tmpi = idx2[ii];
						idx2[ii] = idx2[idxcount-1];
						idx2[idxcount-1] = tmpi;
						idxcount--;
 					}
  					else {
  						q[k]   = tmp;
  						idx[k] = j;
  					}
				}
				
                //If we have only one particle left to choose from then just select it
				if(idxcount<=1){
					for(i=0;i<n;i++){
						xs[n*(t*M+j)+i] = xt[k*n+i];
					}
					numsuccessful++;
				} else if(tmp>1e-16){
                    //Otherwise choose the selected particle with the calculated likelihood
					if (((double) rand())/((double)RAND_MAX) <= tmp){
						for(i=0;i<n;i++){
							xs[n*(t*M+j)+i] = xt[k*n+i];
						}
						numsuccessful++;
					}
				}
			}
		}
	}
	
    
	exit_stub:
		free(x);
		free(q);
		free(epsil);
		free(idx);
		free(idx2);
}
