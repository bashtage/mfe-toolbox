/* Change this to reflect the apropriate header file */
#include <math.h>
#include "mex.h"
#include "matrix.h"

/*
 * igarch_core.c -
 * This is a helper function and is part of the UCSD_GARCH toolbox
 * You can compile it and should work on any platform.
 *
 * Copyright: Kevin Sheppard
 * kevin.sheppard@economics.ox.ac.uk
 * Revision: 3    Date: 9/1/2005
 */
void igarch_core(double *fepsilon, double *parameters, double *backCast, int p, int q, int m, mwSize T, int igarchType, int constant, double *ht) {
    mwIndex i, j;
    
    /* Initialize the final parameter */
    double finalParameter = 1;
    for (i = 0; i < p + q -1; i++) {
        finalParameter -= parameters[constant+i];
    }
    
    for (j=0; j<m; j++) {
        /* uses the estimated cov for starting values */
        ht[j]=backCast[0];
    }
    
    for (i=m; i<T; i++) {
        if (constant) {
            ht[i]=parameters[0];
        }
        for (j=0; j<p; j++) {
            ht[i] += parameters[j+constant] * fepsilon[i-1-j];
        }
        for (j=0; j<q-1; j++) {
            ht[i] += parameters[j+p+constant] * ht[i-1-j];
        }
        ht[i] += finalParameter*ht[i-q];
    }
    if (igarchType==1) {
        for (i=m; i<T; i++) {
            ht[i]=  ht[i] * ht[i];
        }
    }
}


/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /*                    1         2       3      4 5 6 7    8          9
     * ht=igarch_core(fepsilon,parameters,backCast,p,q,m,T,igarchType,constant)
     */
    double *fepsilon, *parameters, *ht, *backCast;
    int p, q, m, igarchType, constant;
    mwSize T;
    /*  Check for proper number of arguments. */
    if(nrhs!=9)
        mexErrMsgTxt("Nine inputs required.");
    if(nlhs>1)
        mexErrMsgTxt("One output only.");
    
    /*  Get the scalar inputs */
    p			= (int)mxGetScalar(prhs[3]);
    q			= (int)mxGetScalar(prhs[4]);
    m			= (int)mxGetScalar(prhs[5]);
    T			= (mwSize)mxGetScalar(prhs[6]);
    igarchType	= (int)mxGetScalar(prhs[7]);
    constant    = (int)mxGetScalar(prhs[8]);
    
    /*  Create a pointer to the input matrices . */
    fepsilon       = mxGetPr(prhs[0]);
    parameters     = mxGetPr(prhs[1]);
    backCast       = mxGetPr(prhs[2]);
    
    /*  Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix(T, 1, mxREAL);
    
    /*  Create a C pointer to a copy of the output matrix. */
    ht = mxGetPr(plhs[0]);
    
    /*  Call the C subroutine. */
    igarch_core(fepsilon, parameters, backCast, p, q, m, T, igarchType, constant, ht);
}