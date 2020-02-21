/* Change this to reflect the apropriate header file */
#include <math.h>
#include <limits.h>
#include "mex.h"
#include "matrix.h"


/*
 * composite_likeliood.c -
 * This is a helper function and is part of the MFE toolbox
 * You should be able to compile it on any platform.
 *
 * Author: Kevin Sheppard
 * kevin.sheppard@economics.ox.ac.uk
 * Revision: 1    Date: 4/17/2012
 */






void composite_likelihood_core(double *S, double *X, double *indices, size_t q, size_t m, size_t n, double *ll) {
    /* q = size(indices,1);
    // [m,n] = size(data);
    // likConst = 3.67575413281869;
    // if m==n
    //     for k=1:q
    //         i = indices(k,1);
    //         j = indices(k,1);
    //         s11 = S(i,i);
    //         s12 = S(i,j);
    //         s22 = S(j,j);
    //         det = s11*s22-s12*s12;
    //         x11 = data(i,i);
    //         x12 = data(i,j);
    //         x22 = data(j,j);
    //         ll = 0.5*(likConst + log(det) + (s22*x11 - 2*s12*x12 + s11*x22)/det)/m;
    //     end
    // else
    //     for k=1:q
    //         i = indices(k,1);
    //         j = indices(k,1);
    //         s11 = S(i,i);
    //         s12 = S(i,j);
    //         s22 = S(j,j);
    //         det = s11*s22-s12*s12;
    //         x11 = data(i)*data(i);
    //         x12 = data(i)*data(j);
    //         x22 = data(j)*data(j);
    //         ll = 0.5*(likConst + log(det) + (s22*x11 - 2*s12*x12 + s11*x22)/det)/m;
    //     end
    // end*/
    size_t k, i, j;
    double s11, s12, s22, x11, x12, x22, scale, det;
    double likConst = 3.67575413281869;
    scale = (double)q;
    *ll = 0.0;
    
    if (m==n)
    {
        for (k=0;k<q;k++)
        {
            i = (size_t)indices[k] - 1;
            j = (size_t)indices[k+q] - 1;
            s11 = S[i*m+i];
            s12 = S[i*m+j];
            s22 = S[j*m+j];
            x11 = X[i*m+i];
            x12 = X[i*m+j];
            x22 = X[j*m+j];
            det = s11*s22 - s12*s12;
            *ll += 0.5*(likConst + log(det) + (s22*x11 - 2*s12*x12 + s11*x22)/det)/scale;
        }
    }
    else
    {
        for (k=0;k<q;k++)
        {
            i = (size_t)indices[k] - 1;
            j = (size_t)indices[k+q] - 1;
            s11 = S[i*m+i];
            s12 = S[i*m+j];
            s22 = S[j*m+j];
            x11 = X[i]*X[i];
            x12 = X[i]*X[j];
            x22 = X[j]*X[j];
            det = s11*s22 - s12*s12;
            *ll += 0.5*(likConst + log(det) + (s22*x11 - 2*s12*x12 + s11*x22)/det)/scale;
        }
    }
}


/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {

    double *S, *X, *indices, *ll;
    size_t m, n, q;
    
    /*  Check for proper number of arguments. */
    if (nrhs!=3)
        mexErrMsgTxt("Three inputs required.");
    if (nlhs>1)
        mexErrMsgTxt("One output only.");
    
    /*  Create a pointer to the input matrices . */
    S       = mxGetPr(prhs[0]);
    X       = mxGetPr(prhs[1]);
    indices = mxGetPr(prhs[2]);
    
    /*  Get the dimensions of the matrix regressand to make an output matrix. */
    m  = (size_t)mxGetM(prhs[1]);
    n =  (size_t)mxGetN(prhs[1]);
    q =  (size_t)mxGetM(prhs[2]);
    
    if (n>m)
    {
        size_t temp;
        temp = m;
        m = n;
        n = temp;
    }
    
    /*  Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    /*  Create a C pointer to a copy of the output matrix. */
    ll = mxGetPr(plhs[0]);
    
    /*  Call the C subroutine. */
    composite_likelihood_core(S, X, indices, q, m, n, ll);
}