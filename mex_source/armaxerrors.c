/* Change this to reflect the apropriate header file */
#include <math.h>
#include <limits.h>
#include "mex.h"
#include "matrix.h"


/*
 * armaxerrors.c -
 * This is a helper function and is part of the UCSD_GARCH/MFE toolbox
 * You should be able to compile it on any platform.
 *
 * Author: Kevin Sheppard
 * kevin.sheppard@economics.ox.ac.uk
 * Revision: 3    Date: 10/16/2009
 */

void armaxerror_core(double *parameters, double *p, double *q, int constant, int np, int nq, double *y, double *x, size_t m, size_t T, size_t k, double *sigma, double *e) {
    
            /*
             * e = zeros(T,1);
             * for t=m+1:T
             * e(t) = y(t);
             * if constant
             * e(t) = e(t) - parameters(1);
             * end
             * for i=1:np
             * e(t) = e(t) - parameters(constant+i)*y(t-p(i));
             * end
             * for i=1:k
             * e(t) = e(t) - parameters(constant+np+i)*x(t,i);
             * end
             * for i=1:nq
             * e(t) = e(t) - parameters(constant+np+k+i)*e(t-q(i));
             * end
             * end
             */
    size_t t, offSet;
    int i;
    for (t = 0; t < m; t++)
    {
        e[t] = 0;
    }
    for (t = m; t < T; t++)
    {
        e[t] = y[t];
        if (constant)
        {
            e[t] -= parameters[0];
        }
        for (i = 0; i < np; i++)
        {
               e[t] -= parameters[constant + i] * y[t - (int)p[i]];
        }
        for (i = 0; i < k; i++)
        {
               offSet = T * i;
               e[t] -= parameters[constant + np + i] * x[t + offSet];
        }
        for (i = 0; i < nq; i++)
        {
               e[t] -= parameters[constant + np + k + i] * e[t - (int)q[i]];
        }
        e[t] = e[t];
    }
    for (t = m; t < T; t++)
    {
        e[t] = e[t]/sigma[t];
    }
}


/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /*                0       1 2 3        4 5 6
     * armaxerrors(parameters,p,q,constant,y,x,m)
     */
    double *parameters, *q, *p, *y, *x, *e, *sigma;
    int constant, np, nq;
    size_t T, m, k;
    
    /*  Check for proper number of arguments. */
    if (nrhs!=8)
        mexErrMsgTxt("Eight inputs required.");
    if (nlhs!=1)
        mexErrMsgTxt("One output only.");
    
    /*  Get the scalar inputs */
    constant = (int)mxGetScalar(prhs[3]);
    m        = (size_t)mxGetScalar(prhs[6]);
    
    /*  Create a pointer to the input matrices . */
    parameters = mxGetPr(prhs[0]);
    p          = mxGetPr(prhs[1]);
    q          = mxGetPr(prhs[2]);
    y          = mxGetPr(prhs[4]);
    x          = mxGetPr(prhs[5]);
    sigma      = mxGetPr(prhs[7]);
    
    /*  Get the dimensions of the matrix regressand to make an output matrix. */
    T  = mxGetM(prhs[4]);
    k =  mxGetN(prhs[5]);
    np =  (int)mxGetM(prhs[1]); 
    nq =  (int)mxGetM(prhs[2]); 

    /*  Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix((mwSize)T, 1, mxREAL);
    
    /*  Create a C pointer to a copy of the output matrix. */
    e = mxGetPr(plhs[0]);
    
    /*  Call the C subroutine. */
    armaxerror_core(parameters, p, q, constant, np, nq, y, x, m, T, k, sigma, e);
}