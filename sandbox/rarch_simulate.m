function [data, Ht] = rarch_simulate(T,C,parameters,p,q,v,type)
% Simulation of RARCH(p,q) multivarate volatility model of Noureldin, Shephard and Sheppard
%
% USAGE:
%  [DATA,HT] = rarch_simulate(T,C,PARAMETERS,P,Q,V,TYPE)
%
% INPUTS:
%   T          - Either a scalar containing the length of the series to simulate, or a T by K matrix 
%                  of simulated random variables.  The default is to use standard normals.  Providing 
%                  a T by K matrix allows other distributions to be used.
%   C          - Unconditional covariance of the data
%   PARAMETERS - Vector of parameters governing the dynamics.  The form of the parameters depends on the TYPE.  
%                  'Scalar':
%                  [a(1) ... a(p) b(1) ... b(q)]'  (all scalars)
%                  'CP' :
%                  [diag(A(:,:,1))' ... diag(A(:,:,p))' theta]' (theta is the scalar persistence)
%                  'Diagonal' 
%                  [diag(A(:,:,1))' ... diag(A(:,:,p))' diag(B(:,:,1))' ... diag(B(:,:,p))']'
%   P          - Positive, scalar integer representing the number of symmetric innovations
%   Q          - Non-negative, scalar integer representing the number of conditional covariance lags
%   V          - Number of components allowed to have time-varying conditional covariance
%   TYPE       - String, one of :
%                  'Scalar' (Default) 
%                  'CP' (Common Persistence) 
%                  'Diagonal'
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
%   and in the CP model, B(:,:,j) = theta - sum(A.^2,3);
%
%   Note that then V<K, then only the first V elements of A and B are
%   non-zero, and PARAMETERS should exclude 0 values.

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012

if strcmpi(type,'Scalar')
    type = 1;
elseif strcmpi(type,'CP')
    type = 2;
elseif strcmpi(type,'Diagonal')
    type = 3;
else
    error('TYPE must be ''Scalar'', ''CP'' or  ''Diagonal''.')
end

k = size(C,1);
if isscalar(T)
    e = randn(2*T,k);
else
    e = T;
    e = [e(ceil(rand(T,1)*T),:);e];
end

[C,A,B] = rarch_parameter_transform(parameters,p,q,v,k,C,type,false,false);

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