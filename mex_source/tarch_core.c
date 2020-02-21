/* Change this to reflect the apropriate header file */
#include <math.h>
#include "mex.h"
#include "matrix.h"

/*
* tarch_core.c -
* This is a helper function and is part of the UCSD_GARCH toolbox
* You can compile it and should work on any platform.
*
* Copyright: Kevin Sheppard
* kevin.sheppard@economics.ox.ac.uk
* Revision: 3    Date: 9/1/2005
*/
void tarch_core(double *fdata, double *fIdata, double *parameters, double *back_cast, int p, int o, int q, int m, mwSize T, int tarch_type, double *ht)
{
	mwIndex i, j;

	for (j=0; j<m; j++)
	{
		/* uses the estimated cov for starting values */
		ht[j]=back_cast[0];
	}

	for (i=m; i<T; i++) {
		ht[i]=parameters[0];
		for (j=0; j<p; j++) {
			ht[i] += parameters[j+1] * fdata[i-1-j];
		}
		for (j=0; j<o; j++) {
			ht[i] += parameters[j+p+1] * fIdata[i-1-j];
		}
		for (j=0; j<q; j++) {
			ht[i] += parameters[j+p+o+1] * ht[i-1-j];
		}
 	}
	if (tarch_type==1)
	{
		for (i=m; i<T; i++)
		{
			ht[i]=  ht[i] * ht[i];
		}
	}
}


/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	/*ht=tarch_core(fdata,fIdata,parameters,back_cast,p,o,q,m,T,tarch_type);*/
	double *fdata, *fIdata, *parameters, *ht, *back_cast;
	int p, o, q, m, tarch_type;
    mwSize T;
	/*  Check for proper number of arguments. */
	if(nrhs!=10)
		mexErrMsgTxt("Ten inputs required.");
	if(nlhs>1)
		mexErrMsgTxt("One output only.");

	/*  Get the scalar inputs */
	p			= (int)mxGetScalar(prhs[4]);
	o			= (int)mxGetScalar(prhs[5]);
	q			= (int)mxGetScalar(prhs[6]);
	m			= (int)mxGetScalar(prhs[7]);
	T			= (mwSize)mxGetScalar(prhs[8]);
	tarch_type	= (int)mxGetScalar(prhs[9]);

	/*  Create a pointer to the input matrices . */
	fdata          = mxGetPr(prhs[0]);
	fIdata         = mxGetPr(prhs[1]);
	parameters     = mxGetPr(prhs[2]);
	back_cast      = mxGetPr(prhs[3]);


	/*  Set the output pointer to the output matrix. */
	plhs[0] = mxCreateDoubleMatrix(T,1, mxREAL);

	/*  Create a C pointer to a copy of the output matrix. */
	ht = mxGetPr(plhs[0]);

	/*  Call the C subroutine. */
	tarch_core(fdata, fIdata, parameters, back_cast, p, o, q, m, T, tarch_type, ht);

}