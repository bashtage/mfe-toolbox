function [O,A,B] = heavy_parameter_transform(parameters,p,q,K)
% Parameter transformation for the HEAVY volatility model of Shephard and Sheppard.  
%
% USAGE:
%  [DATA,HT] = heavy_transform(PARAMETERS,P,Q,K)
%
% INPUTS:
%   PARAMETERS - A vector with K+sum(sum(P))+sum(sum(Q)) elements. See COMMENTS.
%   P          - A K by K matrix containing the lag length of model innovations.  Position (i,j)
%                  indicates the number of lags of series j in the model for series i
%   Q          - A K by K matrix containing the lag length of conditional variances.  Position (i,j)
%                  indicates the number of lags of series j in the model for series i
%   K          - Number of series to simulate
%
% OUTPUTS:
%   O          - A K by 1 matrix of intercepts
%   A          - A K by K  by max(max(P)) array of innovation parameters
%   B          - A K by K  by max(max(Q)) array of smoothing parameters
%
% COMMENTS:
%   Dynamics are given by:
%     h(t,:)' = O + A(:,:,1)*f(data(t-1,:))' + ... + A(:,:,maxP)*f(data(t-maxP,:))' + ...
%                 + B(:,:,1)*h(t-1,:)' + ... + B(:,:,maxQ)*h(t-maxQ,:)'
%
%   PARAMETERS are ordered:
%   [O' A(1,1,1:p(1,1)) A(1,2,1:p(1,2)) ... A(1,K,1:p(1,K)) A(2,1,1:p(2,1)) ... A(2,K,1:p(2,K)) 
%       ... A(K,1,1:p(K,1)) ... A(K,K,1:p(K,K)) ... B(1,1,1:q(1,1)) ... B(1,K,1:q(1,K)) 
%       ... B(K,1,1:q(K,1)) ... B(K,K,1:q(K,K)) ]
%
% EXAMPLES:
%
% See also HEAVY

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/2/2012


O = zeros(K,1);
O(:) = parameters(1:K);
offset = K;
A = zeros(K,K,max(max(p)));
B = zeros(K,K,max(max(q)));
for i=1:K
    for j=1:K
        A(i,j,1:p(i,j)) = parameters(offset+(1:p(i,j)));
        offset = offset + p(i,j);
    end
end
for i=1:K
    for j=1:K
        B(i,j,1:q(i,j)) = parameters(offset+(1:q(i,j)));
        offset = offset + q(i,j);
    end
end
