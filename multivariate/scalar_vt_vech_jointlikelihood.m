function [ll,lls,ht]=scalar_vt_vech_jointlikelihood(parameters,data,p,q,backCast,isJoint)
% Joint (intercept and dymanics parameters) log likelihood for SCALAR_VT_VECH(P,Q) estimation
%
% USAGE:
%   [LL, LLS, HT] = scalar_vt_vech_jointlikelihood(PARAMETERS, DATAAUG, P, Q, K, K2, T)
%
% INPUTS:
%   PARAMETERS    - A vector of vech GARCH process parameters: [vech(C)' alpha beta]'
%   DATAAUG       - Augmented (by m back cast values) matrix of mean zero residuals
%   P             - Positive, scalar integer representing the number of lags of the innovation process
%   Q             - Non-negative scalar integer representing the number of lags of conditional covariance
%   T             - Length of the original data
%
% OUTPUTS:
%   LL             - Minus 1 times the log likelihood
%   LLS            - Time series of log likelihoods (Also multiplied by -1)
%   HT             - Time series of conditional covariances
%
% COMMENTS:
%   See also SCALAR_VT_VECH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

k=size(data,2);
k2=k*(k+1)/2;

C=parameters(1:k2);
C=ivech(C);

%Parse the parameters
alpha=parameters(k2+1:k2+p);
beta=parameters(k2+p+1:k2+p+q);

%Compute the constant
const=(1-sum(parameters))*C;

%Compute the back cast
m=max(p,q);

%Initialize the covariance
Ht=repmat(C,[1 1 t]);

%Initialize the log likelihood
ll=0;
lls=zeros(t+m,1);

%Compute the likelihood constant.
likconst=k*log(2*pi);

%Perform the recursion
for i=m+1:t+m;
    Ht(:,:,i)=const;
    for j=1:p
        Ht(:,:,i)=Ht(:,:,i)+alpha(j)*(data(i-j,:))'*(data(i-j,:));
    end
    for j=1:q
        Ht(:,:,i)=Ht(:,:,i)+beta(j)*Ht(:,:,i-j);
    end
    %Replace these lines to make it work better with poorle conditioned data
    %likelihoods(i)=likconst+(log(det(Ht(:,:,i)))+data(i,:)*Ht(:,:,i)^(-1)*data(i,:)');

    %This is a trick to ensure higher numerical precision
    Q=sqrt(diag(Ht(:,:,i)));
    R=Ht(:,:,i)./(Q*Q');
    stdresid=data(i,:)./Q';
    lls(i)=likconst+2*sum(log(Q))+log(det(R))+stdresid*R^(-1)*stdresid';
    ll=ll+lls(i);
end

%Normalize the log likelihood
ll=0.5*ll;

%If more than one argument out is requested, perform trunctaions and
%normalizations as required
if nargout>1
    lls=0.5*lls(m+1:t+m);
    if nargout>2
        Ht=Ht(:,:,m+1:t+m);
    end
end
