function [ll,lls,Ht]=scalar_vt_vech_likelihood(parameters,data,dataAsym,p,o,q,C,Casym,kappa,backCast,backCastAsym,isJoint,useComposite,estimFlag)
% Log likelihood for SCALAR_VT_VECH(P,O,Q) estimation
%
% USAGE:
%   [LL, LLS, HT] = scalar_vt_vech_likelihood(PARAMETERS,DATA,DATAASYM,P,O,Q,C,CASYM,KAPPA,BACKCAST,BACKCASTASYM,ISJOINT,USECOMPOSITE,ESTIMFLAG)
%
% INPUTS:
%   PARAMETERS   - A vector of vech GARCH process parameters: [alpha beta]' or
%                    [vech(C)' alpha beta]' if ISJOINT
%   DATA         - K by K by T matrix of covariance innovations
%   DATAASYM     - K by K by T matrix of asymmetric covariance innovations
%   P            - Positive, scalar integer representing the number of lags of the innovation process
%   O            - Non-negative scalar integer representing the number of lags of asymmetric process
%   Q            - Non-negative scalar integer representing the number of lags of conditional covariance
%   C            - The unconditional covariance of the data (cov(data)
%   CASYM        - The unconditional expectation of the asymmetric covariance
%   BACKCAST     - Back cast value for starting the recursion
%   BACKCASTASYM - Back cast value (asymetric terms) for starting the recursion
%   ISJOINT      - Logical indicating whether the likelihood is the joint across all parameters
%                    (intercept and dynamics) or is just for the dynamics
%   USECOMPOSITE - Indicates whether to use composite likelihood:
%                    0 - Use standard likelihood (not composite)
%                    1 - Use diagonal composite likelihood
%                    2 - Use full composite likelihood
%   ESTIMFLAG    - [OPTIONAL] Flag (0 or 1) to indicate if the function
%                    is being used in estimation.  If it is 1, then the parameters are transformed
%                    from unconstrained values to constrained by standard garch model constraints
%
% OUTPUTS:
%   LL           - Minus 1 times the log likelihood
%   LLS          - Time series of log likelihoods (Also multiplied by -1)
%   HT           - Time series of conditional covariances
%
% COMMENTS:
%   See also SCALAR_VT_VECH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 4    Date: 10/28/2009

[k,~,T]=size(data);
%If nargin~=7, then we don't need to transform the parameters
base = 0;
if estimFlag
    parameters=scalar_vt_vech_itransform(parameters,p,o,q,kappa);
    %Parse the parameters
elseif isJoint
    k2 = k*(k+1)/2;
    base = k2;
    C = parameters(1:k2);
    C = ivech(C);
    if o>0
        Casym = parameters(k2+1:2*k2);
        Casym = ivech(Casym);
        base = base + k2;
    end
end
alpha=parameters(base+(1:p));
gamma=parameters(base+(p+1:p+o));
beta=parameters(base+(p+o+1:p+o+q));
%Compute the constant
const=(1-sum(alpha)-sum(beta))*C;
if o>0
    const = const - sum(gamma)*Casym;
end

%Initialize the covariance
Ht=repmat(backCast,[1 1 T]);

%Initialize the log likelihood
lls=zeros(T,1);

% Indices or constant, as needed
if useComposite == 0
    likconst = k*log(2*pi);
elseif useComposite == 1
    indices = [(1:k-1)' (2:k)'];
elseif useComposite == 2
    [i,j] = meshgrid(1:k);
    indices = [i(~triu(true(k))) j(~triu(true(k)))];
end
%Perform the recursion
for t=1:T;
    Ht(:,:,t)=const;
    for j=1:p
        if (t-j)<1
            Ht(:,:,t)=Ht(:,:,t)+alpha(j)*backCast;
        else
            Ht(:,:,t)=Ht(:,:,t)+alpha(j)*data(:,:,t-j);
        end
    end
    for j=1:o
        if (t-j)<1
            Ht(:,:,t)=Ht(:,:,t)+gamma(j)*backCastAsym;
        else
            Ht(:,:,t)=Ht(:,:,t)+gamma(j)*dataAsym(:,:,t-j);
        end
    end
    for j=1:q
        if (t-j)<1
            Ht(:,:,t)=Ht(:,:,t)+beta(j)*backCast;
        else
            Ht(:,:,t)=Ht(:,:,t)+beta(j)*Ht(:,:,t-j);
        end
    end
    % Replace these lines to make it work better with poorly conditioned covairance
    % likelihoods(i)=likconst+(log(det(Ht(:,:,t)))+data(i,:)*Ht(:,:,t)^(-1)*data(i,:)');
    % This is a trick to ensure better numerical stability
    if useComposite==0
        Q=sqrt(diag(Ht(:,:,t)));
        R=Ht(:,:,t)./(Q*Q');
        stdresid=data(:,:,t)./(Q*Q');
        lls(t)=0.5*(likconst+2*sum(log(Q))+log(det(R))+trace(R^(-1)*stdresid));
    elseif useComposite
        lls(t) = composite_likelihood(Ht(:,:,t),data(:,:,t),indices);
    end
end
ll = sum(lls);