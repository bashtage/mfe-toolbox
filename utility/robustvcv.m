function [VCV,A,B,scores,hess,gross_scores]=robustvcv(fun,theta,nw,varargin)
% Compute Robust Variance Covariance matrix numerically, including
% Newey-West style score covariance using 2-sided derivatives
%
% USAGE:
%     [VCV,A,B,SCORES,HESS,GROSS_SCORES]=robustvcv(FUN,THETA,NW,VARARGIN)
%
% INPUTS:
%     FUN           - Function name ('fun') or function handle (@fun) which will
%                       return the sum of the log-likelihood (scalar) as the 1st output and the individual
%                       log likelihoods (T by 1 vector) as the second output.
%     THETA         - Parameter estimates at the optimum, usually from fmin*
%     NW            - Number of lags to consider in Newey-West covariance.
%                       Normally set to 0
%     VARARGIN      - Other inputs to the log-likelihood function, such as data
%
% OUTPUTS:
%     VCV           - Estimated robust covariance matrix (see White 1994)
%     A             - A portion of robust covariance
%     B             - B portion of robust covariance
%     SCORES        - T x num_parameters matrix of scores
%     HESS          - Estimated Hessian (Expectation of second derivative)
%     GROSS_SCORES  - Numerical scores (1 by num_parameters) of the objective function, usually for diagnostics
%
% COMMENTS:
%     This function simplifies calculating sandwich covariance estimators for (Q)MLE estimation


% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 9/1/2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(theta,1)<size(theta,2)
    theta=theta';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


k=length(theta);
h=max(abs(theta*eps^(1/3)),1e-8);
h=diag(h);

[~,like]=feval(fun,theta,varargin{:});

t=length(like);

LLFp=zeros(k,1);
LLFm=zeros(k,1);
likep=zeros(t,k);
likem=zeros(t,k);
for i=1:k
    thetaph=theta+h(:,i);
    [LLFp(i),likep(:,i)]=feval(fun,thetaph,varargin{:});
    thetamh=theta-h(:,i);
    [LLFm(i),likem(:,i)]=feval(fun,thetamh,varargin{:});
end

scores=zeros(t,k);
gross_scores=zeros(k,1);
h=diag(h);
for i=1:k
    scores(:,i)=(likep(:,i)-likem(:,i))./(2*h(i));
    gross_scores(i)=(LLFp(i)-LLFm(i))./(2*h(i));
end

hess=hessian_2sided(fun,theta,varargin{:});
A=hess/t;
hess=A;
Ainv=A^(-1);
if nw==0
    % VCV=A^(-1)*B*A^(-1)/t;
    B=cov(scores);
    VCV=(Ainv*B*Ainv)/t;
else
    B=covnw(scores,nw);
    VCV=(Ainv*B*Ainv)/t;
end