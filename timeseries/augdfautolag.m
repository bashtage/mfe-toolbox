function [adfstat,pval,critval,resid,lags,ICs]=augdfautolag(y,p,maxlags,IC)
% Dickey-Fuller and Augmented Dickey Fuller with automatic lag selection
%
% USAGE:
%  [ADFSTAT,PVAL,CRITVAL] = augdfautolag(Y,P,LAGS,IC)
%  [ADFSTAT,PVAL,CRITVAL,RESID,LAGS] = augdfautolag(Y,P,LAGS,IC)
%
% INPUTS:
%  Y         - A T by 1 vector of data
%  P         - Order of the polynomial of include in the ADF regression:
%                0 : No deterministic terms
%                1 : Constant
%                2 : Time Trend
%                3 : Constant, DGP assumed to have a time trend
%  MAXLAGS   - The maximum number of lags to include in the ADF test
%  IC        - [OPTIONAL] String, either 'AIC' (default) or 'BIC' to choose the criteria to select
%                the model
%
% OUTPUTS:
%  ADFSTAT   - Dickey-Fuller statistic
%  PVAL      - Probability the series is a unit root
%  CRITVALS  - A 6 by 1 vector with the [.01 .05 .1 .9 .95 .99] values from the DF distribution
%  LAGS      - The selected number of lags
%  IC        - The value at all lags of the selected IC
%
% COMMENTS:
%
% See also AUGDF

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3.0.1    Date: 1/1/2007



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3 || nargin>4
    error('3 or 4 inputs required.')
end
T=length(y);
if T<=(maxlags+1)
    error('Length of data must be larger than LAGS')
end
if size(y,1)~=T,
    y=y';
end
if size(y,2)~=1
    error('Y must be a column vector')
end
if ~isscalar(maxlags) && maxlags>0 && floor(maxlags)==maxlags
    error('LAGS must be a positive integer')
end
if ~isscalar(p)
    error('P must be a scalar integer in {0, 1, 2, 3}')
end
if ~ismember(p,[0 1 2 3])
    error('P must be a scalar integer in {0, 1, 2, 3}')
end
if nargin==3
    IC='AIC';
else
    if ~ischar(IC)
        error('IC must be a string, either ''AIC'' or ''BIC''');
    elseif ~ismember(IC,{'AIC','BIC'})
        error('IC must be a string, either ''AIC'' or ''BIC''');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Setup common to all problems
%y=y-y(1);
ydiff=diff(y);
[ydiffcurr, ydifflags]=newlagmatrix(ydiff,maxlags); %#ok<ASGLU>
T=length(y);
Y=y(maxlags+2:T);
tau=length(Y);
%
switch p
    case 0
        % Case 1
        i=0;
        X=y(maxlags+1:T-1);
        rho = X\Y;
        % Compute the errors
        e= Y-X*rho;
        s2(i+1)=e'*e/tau;
        K(i+1)=size(X,2);
        
        % Loop
        for i=1:maxlags
            X=[y(maxlags+1:T-1) ydifflags(:,1:i)];
            rho = X\Y;
            % Compute the errors
            e= Y-X*rho;
            s2(i+1)=e'*e/tau;
            K(i+1)=size(X,2);
        end
        
    case {1,3}
        %Case 2
        i=0;
        X=[ones(size(Y)) y(maxlags+1:T-1)];
        rho = X\Y;
        % Compute the errors
        e= Y-X*rho;
        s2(i+1)=e'*e/tau;
        K(i+1)=size(X,2);
        
        % Loop
        for i=1:maxlags
            X=[ones(size(Y)) y(maxlags+1:T-1) ydifflags(:,1:i)];
            rho = X\Y;
            % Compute the errors
            e= Y-X*rho;
            s2(i+1)=e'*e/tau;
            K(i+1)=size(X,2);
        end
        
        
    case 2
        %Case 4
        i=0;
        X=[ones(size(Y)) y(maxlags+1:T-1) (1:tau)'];
        rho = X\Y;
        % Compute the errors
        e= Y-X*rho;
        s2(i+1)=e'*e/tau;
        K(i+1)=size(X,2);
        
        % Loop
        for i=1:maxlags
            X=[ones(size(Y)) y(maxlags+1:T-1) (1:tau)' ydifflags(:,1:i)];
            rho = X\Y;
            % Compute the errors
            e= Y-X*rho;
            s2(i+1)=e'*e/tau;
            K(i+1)=size(X,2);
        end
end


if strcmp(IC,'AIC')
    ICs=log(s2) + 2*K/tau;
else
    ICs=log(s2) + K*log(tau)/tau;
end
[~,lags]=min(ICs);
lags=lags-1;
[adfstat,pval,critval,resid]=augdf(y,p,lags);