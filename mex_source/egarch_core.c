/* Change this to reflect the apropriate header file */
#include <math.h>
#include <limits.h>
#include "mex.h"
#include "matrix.h"

/*
* egarch_core.c -
* This is a helper function and is part of the UCSD_GARCH toolbox
* You can compile it and should work on any platform.
*
* Author: Kevin Sheppard
* kevin.sheppard@economics.ox.ac.uk
* Revision: 3    Date: 3/1/2005
*/
void egarch_core(double *data, double *parameters, double *back_cast, double *upper, int p, int o, int q, int m, mwSize T, double *ht, double *logHt, double *stdData, double *absStdData)
{
    mwIndex i;
	int j;
	double exp_back_cast, vol, subconst, eLB, hLB;
	
	exp_back_cast=exp(back_cast[0]);
	subconst = 0.797884560802865;
	for (j=0; j<m; j++) 
	{
		/* uses the estimated cov for starting values*/
		logHt[j]=back_cast[0];
		ht[j]=exp_back_cast;
		vol = sqrt(ht[j]);
		stdData[j] = data[j]/vol;
		absStdData[j] = fabs(stdData[j])-subconst;
	}
    eLB = exp_back_cast/10000;
    hLB = back_cast[0] - log(10000);
	for (i=m; i<T; i++) {
		logHt[i]=parameters[0];
		for (j=0; j<p; j++) {
			logHt[i] += parameters[j+1] * absStdData[i-1-j];
		}
		for (j=0; j<o; j++) {
			logHt[i] += parameters[j+p+1] * stdData[i-1-j];
		}
		for (j=0; j<q; j++) {
			logHt[i] += parameters[j+p+o+1] * logHt[i-1-j];
		}
		/* Compute the stanardized residuals*/
		ht[i]=exp(logHt[i]);
		if (ht[i]<eLB)
		{
			ht[i]=eLB;
			logHt[i] = hLB;
		}
		/*Should have a check that not to big*/
		if (ht[i]>upper[0])
		{
			if(mxIsInf(ht[i]))
			{
				ht[i]  = upper[0];
			}
			else
			{
				ht[i]= upper[0] + logHt[i];
			}
			logHt[i] = log(upper[0]);
		}
		vol = sqrt(ht[i]);
		stdData[i] = data[i]/vol;
		absStdData[i] = fabs(stdData[i])-subconst;
 	}
}


/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	/*ht=egarch_core(data,parameters,back_cast,upper,p,o,q,m,T);*/
	double *data, *parameters, *ht, *back_cast, *upper, *logHt, *stdData, *absStdData;
	int p, o, q, m;
    mwSize T;
	mxArray *pTemp[3];

	/*  Check for proper number of arguments. */
	if(nrhs!=9)
		mexErrMsgTxt("Nine inputs required.");
	if(nlhs>1)
		mexErrMsgTxt("One output only.");

	/*  Get the scalar inputs */
	p			= (int)mxGetScalar(prhs[4]);
	o			= (int)mxGetScalar(prhs[5]);
	q			= (int)mxGetScalar(prhs[6]);
	m			= (int)mxGetScalar(prhs[7]);
	T			= (mwSize)mxGetScalar(prhs[8]);

	/*  Create a pointer to the input matrices . */
	data           = mxGetPr(prhs[0]);
	parameters     = mxGetPr(prhs[1]);
	back_cast      = mxGetPr(prhs[2]);
	upper		   = mxGetPr(prhs[3]);

	
	/*  Set the output pointer to the output matrix. */
	plhs[0] = mxCreateDoubleMatrix(T,1, mxREAL);
	pTemp[0] = mxCreateDoubleMatrix(T,1, mxREAL);
	pTemp[1] = mxCreateDoubleMatrix(T,1, mxREAL);
	pTemp[2] = mxCreateDoubleMatrix(T,1, mxREAL);

	/*  Create a C pointer to a copy of the output matrix. */
	ht = mxGetPr(plhs[0]);
	logHt = mxGetPr(pTemp[0]);
	stdData = mxGetPr(pTemp[1]);
	absStdData = mxGetPr(pTemp[2]);

	/*  Call the C subroutine. */
	egarch_core(data, parameters, back_cast, upper, p, o, q, m, T, ht, logHt, stdData, absStdData);
	mxDestroyArray(pTemp[0]);
	mxDestroyArray(pTemp[1]);
	mxDestroyArray(pTemp[2]);
}