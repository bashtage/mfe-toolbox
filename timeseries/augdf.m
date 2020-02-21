function [adfstat,pval,critval,resid]=augdf(y,p,lags)
% Dickey-Fuller and Augmented Dickey Fuller testing
%
% USAGE:
%  [ADFSTAT,PVAL,CRITVAL] = augdf(Y,P,LAGS)
%  [ADFSTAT,PVAL,CRITVAL,RESID] = augdf(Y,P,LAGS)
% 
% INPUTS:
%  Y         - A T by 1 vector of data
%  P         - Order of the polynomial of include in the ADF regression:
%                0 : No deterministic terms
%                1 : Constant
%                2 : Time Trend
%                3 : Constant, DGP assumed to have a time trend
%  LAGS      - The number of lags to include in the ADF test (0 for DF test)
% 
% OUTPUTS:
%  ADFSTAT   - Dickey-Fuller statistic
%  PVAL      - Probability the series is a unit root
%  CRITVALS  - A 6 by 1 vector with the [.01 .05 .1 .9 .95 .99] values from the DF distribution
%  RESID     - Residual (adjusted for lags) from the ADF regression
% 
% COMMENTS:
%  
% See also AUGDFAUTOLAG


% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3.0.1    Date: 1/1/2007



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=3
    error('3 inputs required.')
end
T=length(y);
if T<=(lags+1)
    error('Length of data must be larger than LAGS')
end
if size(y,1)~=T,
    y=y';
end
if size(y,2)~=1
    error('Y must be a column vector')
end
if ~isscalar(lags) && lags>0 && floor(lags)==lags
    error('LAGS must be a positive integer')
end
if ~isscalar(p)
    error('P must be a scalar integer in {0, 1, 2, 3}')
end
if ~ismember(p,[0 1 2 3])
    error('P must be a scalar integer in {0, 1, 2, 3}')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Setup common to all problems
%y=y-y(1);
ydiff=diff(y);
[ydiffcurr, ydifflags]=newlagmatrix(ydiff,lags);
T=length(y);
Y=y(lags+2:T);
tau=length(Y);
%
switch p
    case 0
        % Case 1
        X=[y(lags+1:T-1) ydifflags];
        rho = X\ydiffcurr;
        % Compute the errors
        e= ydiffcurr-X*rho;

        % Compute the covariance matrix of the estimated parameters
        s2 = e'*e/(tau-size(X,2));
        Uinv=inv(diag([T T^(0.5)*ones(1,lags)]));
        sel=[1 zeros(1,lags)];
        sigp=s2*sel*(Uinv*(X'*X)*Uinv)^(-1)*sel';

        % Compute the ADF
        adfstat = tau*rho(1)/sqrt(sigp);
        % Look up the pval and critical values
        [pval,critval]=augdfcv(adfstat,p,tau);
    case 1
        %Case 2
        X=[ones(size(Y)) y(lags+1:T-1) ydifflags];
        rho = X\ydiffcurr;
        % Compute the errors
        e= ydiffcurr-X*rho;

        % Compute the covariance matrix of the estimated parameters
        s2 = e'*e/(tau-size(X,2));
        Uinv=inv(diag([T^(0.5) T T^(0.5)*ones(1,lags)]));
        sel=[0 1 zeros(1,lags)];
        sigp=s2*sel*(Uinv*(X'*X)*Uinv)^(-1)*sel';

        % Compute the ADF
        adfstat = tau*rho(2)/sqrt(sigp);
        % Look up the pval and critical values
        [pval,critval]=augdfcv(adfstat,p,tau);

    case 2
        %Case 4
        X=[ones(size(Y)) y(lags+1:T-1) (1:tau)' ydifflags];
        rho = X\Y;
        % Compute the errors
        e= Y-X*rho;

        % Compute the covariance matrix of the estimated parameters
        s2 = e'*e/(tau-size(X,2));
        X(:,2)=X(:,2)-rho(1)*(1:tau)';
        Uinv=inv(diag([tau^(0.5) tau tau^(1.5) tau^(0.5)*ones(1,lags)]));
        sel=[0 1 0 zeros(1,lags)];
        sigp=s2*sel*(Uinv*(X'*X)*Uinv)^(-1)*sel';

        % Compute the ADF
        adfstat = tau*(rho(2)-1)/sqrt(sigp);
        % Look up the pval and critical values
        [pval,critval]=augdfcv(adfstat,p,tau);

    case 3
        %Case 3
        X=[ones(size(Y)) y(lags+1:T-1) ydifflags];
        rho = X\Y;
        % Compute the errors
        e= Y-X*rho;

        % Compute the covariance matrix of the estimated parameters
        s2 = e'*e/(tau-size(X,2));
        Uinv=inv(diag([tau^(0.5) tau^(1.5)  tau^(0.5)*ones(1,lags)]));
        sel=[0 1 zeros(1,lags)];
        sigp=s2*sel*(Uinv*(X'*X)*Uinv)^(-1)*sel';

        % Compute the ADF
        adfstat = tau^(1.5)*(rho(2)-1)/sqrt(sigp);
        % Look up the pval and critical values
        critval=norminv([.01 .05 .1 .9 .95 .99]');
        pval = normcdf(adfstat);
end
resid=e;