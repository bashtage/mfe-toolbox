function [errors,E]=armaxfilter_core(regressand,regressors, arxparameters, maparameters, q, K, tau, maxq, maxe)
% M-file version of the inner loop in conputing the forward recursion for
% MA, MAX, ARMA and ARMAX estimation
% 
% USAGE:
%     [E,E2]=armaxfilter_core(REGRESSAND,REGRESSORS,ARXPARAMETERS,MAPARAMETERS,Q,TAU,MAXQ,MAXE)
% 
% INPUTS:
%     See armaxlikelihood
% 
% OUTPUTS:
%     See armaxlikelihood
% 
% COMMENTS:
%     Helper function part of UCSD_GARCH toolbox. Used if you do not use the MEX file.
%     You should use the MEX file.
% 
% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 1/1/2007

errors=zeros(tau,1);
nq = size(q,1);
for t=(1+maxq):tau;
    errors(t)=regressand(t);
    for i=1:K
        errors(t) = errors(t) - arxparameters(i)*regressors(t,i);
    end
    for i=1:nq
        errors(t)=  errors(t) - maparameters(i)*errors(t-q(i));
    end
    if abs(errors(t))>(100000*maxe)
        errors(t)=maxe*sign(errors(t));
    end
end
E=errors(maxq+1:tau)'*errors(maxq+1:tau);
