function [data, Ht] = bekk_simulate(T,k,parameters,p,o,q,type)
% Simulation of symmetric and asymmetric BEKK(p,o,q) multivariate volatility models
%
% USAGE:
%  [DATA,HT] = bekk_simulate(T,K,PARAMETERS,P,O,Q,TYPE)
%
% INPUTS:
%   T          - Either a scalar containing the length of the series to simulate, or a T by K matrix 
%                  of simulated random variables.  The default is to use standard normals.  Providing 
%                  a T by K matrix allows other distributions to be used.
%   PARAMETERS - Vector of parameters.  The form of the parameters depends on the TYPE.  
%                  'Scalar':
%                  [CC' a(1) ... a(p) g(1) ... g(o) b(1) ... b(q)]'  (all scalars)
%                  'Diagonal' 
%                  [CC' diag(A(:,:,1))' ... diag(A(:,:,p))' diag(G(:,:,1))' ... diag(G(:,:,o))' diag(B(:,:,1))' ... diag(B(:,:,p))']'
%                  'Full' 
%                  [CC' f(A(:,:,1)) ... f(A(:,:,p)) f(G(:,:,1)) ... f(G(:,:,o)) f(B(:,:,1)) ... f(B(:,:,q))]'
%                  where CC = chol2vec(C')' and f(M) = M(:)'
%   P          - Positive, scalar integer representing the number of symmetric innovations
%   O          - Non-negative, scalar integer representing the number of asymmetric innovations
%   Q          - Non-negative, scalar integer representing the number of conditional covariance lags
%   TYPE       - String, one of :
%                  'Scalar' (Default) 
%                  'Diagonal'
%                  'Full'
%
% OUTPUTS:
%   DATA   - A T by K matrix of simulated data
%   HT     - A [K K T] dimension matrix of conditional covariances
%
% COMMENTS:
%   The dynamics of a BEKK are given by 
%   
%   H(:,:,t) = C*C' +
%       A(:,:,1)'*OP(:,:,t-1)*A(:,:,1) + ... + A(:,:,p)'*OP(:,:,t-1)*A(:,:,p) +
%       G(:,:,1)'*OPA(:,:,t-1)*G(:,:,1) + ... + G(:,:,o)'*OPA(:,:,t-1)*G(:,:,o) +
%       B(:,:,1)'*G(:,:,t-1)*B(:,:,1) + ... + B(:,:,q)'*OP(:,:,t-1)*B(:,:,q)
%
%   where in the scalar model A(:,:,i) = a(i)*eye(K) (similarly for G and B).
%
%  EXAMPLES:
%    % Scalar with A.^2=.05, G.^2=.1 and B.^2=.88
%    CCp = [1 .5;.5 4];
%    parameters = [chol2vec(chol(CCp)');sqrt([.05,.10,.88])']
%    [data,Ht] = bekk_simulate(1000,2,parameters,1,1,1,'Scalar')
%    % Diagonal 
%    parameters = [chol2vec(chol(CCp)');sqrt([.05 .07 .93 .88])']
%    [data,Ht] = bekk_simulate(1000,2,parameters,1,0,1,'Diagonal')
%
% See also BEKK


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
elseif strcmpi(type,'Diagonal')
    type = 2;
elseif strcmpi(type,'Full')
    type = 3;
else
    error('TYPE must be ''Scalar'', ''Diagonal'' or ''Full''.')
end

if p<=0 || floor(p)~=p
    error('P must be a positive scalar integer.')
end
if o<0 || floor(o)~=o
    error('O must be a positive scalar integer.')
end
if q<0 || floor(q)~=q
    error('Q must be a positive scalar integer.')
end

k2 = k*(k+1)/2;
switch type
    case 1
        count = p+o+q;
    case 2
        count = (p+o+q)*k;
    case 3
        count = (p+o+q)*k*k;
end
count = count + k2;
if length(parameters)~=count
    error('PARAMETERS does not have the expected number of elements.')
end
[C,A,G,B] = bekk_parameter_transform(parameters,p,o,q,k,type);
m = zeros(k*k);
for i=1:p
    m = m + kron(A(:,:,i),A(:,:,i));
end
for i=1:o
    m = m + 0.5*kron(G(:,:,i),G(:,:,i));
end
for i=1:q
    m = m + kron(B(:,:,i),B(:,:,i));
end

if max(eigs(m))>1
    backCast = C/.001;
    warning('MFE:nonstationary','The parameters do not correspond to the stationary region.')
else
    backCast = ((eye(k*k)-m)\eye(k*k))*C(:);
    backCast = reshape(backCast,k,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ht = repmat(eye(k),[1 1 2*T]);
eta = zeros(size(e));
for i=1:(2*T)
    Ht(:,:,i) = C;
    for j=1:p
        if (i-j)<=0
            Ht(:,:,i) = Ht(:,:,i) + A(:,:,j)'*backCast*A(:,:,j);
        else
            Ht(:,:,i) = Ht(:,:,i) + A(:,:,j)'*(e(i-j,:)'*e(i-j,:))*A(:,:,j);
        end
    end
    for j=1:o
        if (i-j)<=0
            Ht(:,:,i) = Ht(:,:,i) + G(:,:,j)'*backCast*G(:,:,j);
        else
            Ht(:,:,i) = Ht(:,:,i) + G(:,:,j)'*(eta(i-j,:)'*eta(i-j,:))*G(:,:,j);
        end
    end    
    for j=1:q
        if (i-j)<=0
            Ht(:,:,i) = Ht(:,:,i) + B(:,:,j)'*backCast*B(:,:,j);
        else
            Ht(:,:,i) = Ht(:,:,i) + B(:,:,j)'*Ht(:,:,i-j)*B(:,:,j);
        end
    end
    Ht12 = Ht(:,:,i)^(0.5);
    e(i,:) = e(i,:)*Ht12;
    eta(i,:) = e(i,:).*(e(i,:)<0);
end

data = e(T+1:2*T,:);
Ht = Ht(:,:,T+1:2*T);