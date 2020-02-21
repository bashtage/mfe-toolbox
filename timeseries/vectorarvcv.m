function VCV = vectorarvcv(X,errors,het,uncorr)
% Estimate the variance-covariance matrix for parameters estimated using vectorar under one of these
% conditions 
%   * Conditionally Homoskedastic and Uncorrelated
%   * Conditionally Homoskedastic but Correlated
%   * Heteroskedastic but Conditionally Uncorrelated
%   * Heteroskedastic and Correlated
%
% USAGE:
%   [VCV] = vectorarvar(X,ERRORS,HET,UNCORR)
%
% INPUTS:
%   X      - A T by # parameters matrix of regressors for each Y
%   ERRORS - K by T vector of errors
%   HET    - A scalar integer indicating the type of covariance estimator
%              0 - Homoskedastic 
%              1 - Heteroskedastic [DEFULAT]
%   HET    - A scalar integer indicating the assumed structure of the error covariance matrix
%              0 - Correlated errors  [DEFULAT]
%              1 - Uncorrelated errors
%
% OUTPUTS:
%   VCV     - A K*((# lags) + CONSTANT) by K*((# lags) + CONSTANT) matrix of estimated parameter
%               covariances computed using  HET and UNCORR 
%
% COMMENTS:
%   Helper function for VECTORAR
%
% See also VECTORAR, GRANGERCAUSE
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3.0    Date: 1/1/2007


T = size(X,1);
s2=errors'*errors/T;
Np = size(X,2);
XpXi = ((X'*X)/T)^(-1);
K=size(errors,2);
X2=repmat(X,1,K);
e2=reshape(repmat(errors,Np,1),T,K*Np);
s=X2.*e2;

if het==1 && uncorr==0
    % This is the White Cov Matrix
    Ainv = kron(eye(K),XpXi);
    B=(s'*s)/T;
    VCV=Ainv*B*Ainv/T;
elseif het==1 && uncorr==1
    Ainv = kron(eye(K),XpXi);    
    B=zeros(Np*K);
    for i=1:K
        sel=(i-1)*Np+1:i*Np;
        temp = s(:,sel);
        B(sel,sel)=temp'*temp/T;
    end
    VCV=Ainv*B*Ainv/T;
elseif het==0 && uncorr==0
    % Homoskedastic, correlated is also fairly easy
    VCV=kron(s2,XpXi)/T;
elseif het==0 && uncorr==1
    % Hom uncorrelated
    VCV=kron(diag(diag(s2)),XpXi)/T;
end
% White uncorr is also easy, since it comes from B
