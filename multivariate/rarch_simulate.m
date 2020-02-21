function [data, Ht] = rarch_simulate(T,C,parameters,p,q,type)
% Simulation of RARCH(p,q) multivariate volatility model of Noureldin, Shephard and Sheppard
%
% USAGE:
%  [DATA,HT] = rarch_simulate(T,C,PARAMETERS,P,Q,TYPE)
%
% INPUTS:
%   T          - Either a scalar containing the length of the series to simulate, or a T by K matrix 
%                  of simulated random variables.  The default is to use standard normal random 
%                  variables.  Providing a T by K matrix allows other distributions to be used.
%   C          - Unconditional covariance of the data
%   PARAMETERS - Vector of parameters governing the dynamics.  The form of the parameters depends on the TYPE.  
%                  'Scalar':
%                  [a(1) ... a(p) b(1) ... b(q)]'  (all scalars)
%                  'CP' :
%                  [diag(A(:,:,1))' ... diag(A(:,:,p))' theta]' (theta is the scalar persistence)
%                  'Diagonal' 
%                  [diag(A(:,:,1))' ... diag(A(:,:,p))' diag(B(:,:,1))' ... diag(B(:,:,p))']'
%   P          - Positive, scalar integer representing the number of symmetric innovations
%   Q          - Non-negative, scalar integer representing the number of conditional covariance
%                  lags.  When using 'CP' model, 0<=q<=1
%   TYPE       - String, one of :
%                  'Scalar' (Default) 
%                  'CP' (Common Persistence) 
%                  'Diagonal'
%
% OUTPUTS:
%   DATA   - A T by K matrix of simulated data
%   HT     - A [K K T] dimension matrix of conditional covariances
%
% COMMENTS:
%   The dynamics of a RARCH model are identical to that of a BEKK, except
%   that the model evolves in the rotated space.
%   
%   G(:,:,t) = (eye(K) - sum(A.^2,3) - sum(B.^2,3)) +
%       A(:,:,1)*OP(:,:,t-1)*A(:,:,1) + ... A(:,:,p)*OP(:,:,t-1)*A(:,:,p) +
%       B(:,:,1)*G(:,:,t-1)*B(:,:,1) + ... B(:,:,p)*OP(:,:,t-1)*B(:,:,p)
%
%   where in the scalar model A(:,:,i) = a(i)*eye(K), B(:,:,j)=b(j)*eye(K)
%   and in the CP model, B(:,:,j) = theta - sum(A.^2,3).  OP is the outer product of the
%   unconditionally standardized data.
%
% EXAMPLES:
%   % Scalar with A.^2=.05 and B.^2=.93
%   [data,Ht] = rarch_simulate(1000,eye(2)+1,sqrt([.05,.93]),1,1,'Scalar')
%   % Diagonal 
%   [data,Ht] = rarch_simulate(1000,eye(2)+1,sqrt([.05 .07 .93 .88]),1,1,'Diagonal')
%   % Common Persistence, note uses theta (sqrt(A.^2+B.^2) not B)
%   [data,Ht] = rarch_simulate(1000,eye(2)+1,sqrt([.05 .07 .99]),1,1,'CP')
%
% See also RARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = size(C,1);
if isscalar(T)
    e = randn(2*T,k);
else
    e = T;
    if size(e,2)~=k
        error('T must have K columns when providing simulated random numbers.')
    end
    T = size(e,1);
    e = [e(ceil(rand(T,1)*T),:);e];
end

if strcmpi(type,'Scalar')
    type = 1;
elseif strcmpi(type,'CP')
    type = 2;
elseif strcmpi(type,'Diagonal')
    type = 3;
else
    error('TYPE must be ''Scalar'', ''CP'' or  ''Diagonal''.')
end

if type==2 && q>1
    error('Q must be either 0 or 1 for the ''CP'' model.')
end

switch type
    case 1
        count = p+q;
    case 2
        count = p*k+q;
    case 3
        count = (p+q)*k;
end
if length(parameters)~=count
    error('PARAMETERS does not have the expected number of elements.')
end
[C,A,B] = rarch_parameter_transform(parameters,p,q,k,C,type,false);
if max(diag(sum(A.^2,3)+sum(B.^2,3)))>=1
    warning('MFE:nonstationary','The parameters do not correspond to the stationary region.')
end
if type==2 && q==1 && min(B(:))<0
    error('When using ''CP'', the common persistence parameter Theta must satisfy Theta^2>=sum(A.^2,3)')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gt = repmat(eye(k),[1 1 2*T]);
intercept = eye(k) - sum(A.^2,3) - sum(B.^2,3);
backCast = eye(k);
for i=1:(2*T)
    Gt(:,:,i) = intercept;
    for j=1:p
        if (i-j)<=0
            Gt(:,:,i) = Gt(:,:,i) + A(:,:,j)*backCast*A(:,:,j);
        else
            Gt(:,:,i) = Gt(:,:,i) + A(:,:,j)*(e(i-j,:)'*e(i-j,:))*A(:,:,j);
        end
    end
    for j=1:q
        if (i-j)<=0
            Gt(:,:,i) = Gt(:,:,i) + B(:,:,j)*backCast*B(:,:,j);
        else
            Gt(:,:,i) = Gt(:,:,i) + B(:,:,j)*Gt(:,:,i-j)*B(:,:,j);
        end
    end
    Gt12 = Gt(:,:,i)^(0.5);
    e(i,:) = e(i,:)*Gt12;
end

data = zeros(size(e));
Ht = zeros(k,k,T);
C12 = C^(0.5);
for i=1:length(e)
    data(i,:) = e(i,:)*C12;
    Ht(:,:,i) = C12*Gt(:,:,i)*C12;
end

data = data(T+1:2*T,:);
Ht = Ht(:,:,T+1:2*T);