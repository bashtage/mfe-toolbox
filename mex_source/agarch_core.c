/* Change this to reflect the apropriate header file */
#include <math.h>
#include "mex.h"
#include "matrix.h"

/*
 * agarch_core.c -
 * This is a helper function and is part of the UCSD_GARCH toolbox
 * You can compile it and should work on any platform.
 *
 * Copyright: Kevin Sheppard
 * kevin.sheppard@economics.ox.ac.uk
 * Revision: 1    Date: 7/12/2009
 */
void agarch_core(double *data, double *parameters, double *back_cast, int p, int q, int m, mwSize T, int model_type, double *ht, double *shock) {

    mwIndex i, j;
    double gamma, gamma2;
    gamma = parameters[p+1];
    gamma2 = gamma * gamma;

    if (model_type == 1) {
        for (j=0; j<m; j++) {
            /* uses the estimated cov for starting values */
            ht[j] = back_cast[0];
            shock[j] = back_cast[0] + gamma2;
        }

        for (i=m; i<T; i++) {
            ht[i]=parameters[0];
            for (j=0; j<p; j++) {
                ht[i] += parameters[j+1] * shock[i-1-j];
            }
            for (j=0; j<q; j++) {
                ht[i] += parameters[j+p+2] * ht[i-1-j];
            }
            shock[i] = data[i] - gamma;
            shock[i] = shock[i] * shock[i];
        }
    }
    else {
        for (j=0; j<m; j++) {
            /* uses the estimated cov for starting values */
            ht[j] = back_cast[0];
            shock[j] = back_cast[0] * (1 + gamma2);
        }

        for (i=m; i<T; i++) {
            ht[i]=parameters[0];
            for (j=0; j<p; j++) {
                ht[i] += parameters[j+1] * shock[i-1-j];
            }
            for (j=0; j<q; j++) {
                ht[i] += parameters[j+p+2] * ht[i-1-j];
            }
            shock[i] = data[i] - gamma * sqrt(ht[i]);
            shock[i] = shock[i] * shock[i];
        }
    }
}


/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    /*ht=agarch_core(data,parameters,back_cast,p,q,m,T,model_type)*/
    double *data, *parameters, *ht, *back_cast, *shock;
    int p, q, m, model_type;
    mwSize T;
    mxArray *pTemp[1];

    /*  Check for proper number of arguments. */
    if(nrhs!=8)
        mexErrMsgTxt("Eight inputs required.");
    if(nlhs>1)
        mexErrMsgTxt("One output only.");

    /*  Get the scalar inputs */
    p			= (int)mxGetScalar(prhs[3]);
    q			= (int)mxGetScalar(prhs[4]);
    m			= (int)mxGetScalar(prhs[5]);
    T			= (mwSize)mxGetScalar(prhs[6]);
    model_type	= (int)mxGetScalar(prhs[7]);

    /*  Create a pointer to the input matrices . */
    data           = mxGetPr(prhs[0]);
    parameters     = mxGetPr(prhs[1]);
    back_cast      = mxGetPr(prhs[2]);


    /*  Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix(T, 1, mxREAL);
    pTemp[0] = mxCreateDoubleMatrix(T,1, mxREAL);

    /*  Create a C pointer to a copy of the output matrix. */
    ht = mxGetPr(plhs[0]);
    shock = mxGetPr(pTemp[0]);


    /*  Call the C subroutine. */
    agarch_core(data, parameters, back_cast, p, q, m, T, model_type, ht, shock);

    /* Clean up */
    mxDestroyArray(pTemp[0]);
}