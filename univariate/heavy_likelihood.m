function [ll,lls,h] = heavy_likelihood(parameters,data,p,q,backCast,lb,ub)
% Likelihood for HEAVY volatility model of Shephard and Sheppard
%
% USAGE:
%  [LL,LLS,H] = heavy_likelihood(PARAMETERS,DATA,P,Q,BACKCAST,LB,UB)
%
% INPUTS:
%   PARAMETERS - A vector with K+sum(sum(P))+sum(sum(Q)) elements. See COMMENTS.
%   DATA       - A T by K vector of non-negative data.  Returns should be squared before using
%   P          - A K by K matrix containing the lag length of model innovations.  Position (i,j)
%                  indicates the number of lags of series j in the model for series i
%   Q          - A K by K matrix containing the lag length of conditional variances.  Position (i,j)
%                  indicates the number of lags of series j in the model for series i
%   BACKCAST   - A 1 by K matrix of values to use fo rback casting
%   LB         - A 1 by K matrix of volatility lower bounds to use in estimation
%   UB         - A 1 by K matrix of volatility upper bounds to use in estimation
%
% OUTPUTS:
%   LL          - The log likelihood evaluated at the PARAMETERS
%   LLS         - A T by 1 vector of log-likelihoods
%   HT          - A T by K matrix of conditional variances
%
% COMMENTS:
%
% EXAMPLES:
%
% See also HEAVY

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/2/2012

pMax = max(max(p));
qMax = max(max(q));
[K,T] = size(data);
lls = zeros(T,1);
h = zeros(K,T);

[O,A,B] = heavy_parameter_transform(parameters,p,q,K);
O(O<0) = realmin();

likConst = K*log(2*pi);
for t=1:T
    h(:,t)= O;
    for j=1:pMax
        if (t-j)>0
            h(:,t) = h(:,t) + A(:,:,j)*data(:,t-j);
        else
            h(:,t) = h(:,t) + A(:,:,j)*backCast';
        end
    end
    for j=1:qMax
        if (t-j)>0
            h(:,t) = h(:,t) + B(:,:,j)*h(:,t-j);
        else
            h(:,t) = h(:,t) + B(:,:,j)*backCast';
        end
    end
    for j=1:K
        if h(j,t)<lb(j)
            h(j,t) = lb(j) * 1./(1-(h(j,t)-lb(j)));
        elseif h(j,t)>ub(j)
            h(j,t) = ub(j) + log(h(j,t)-ub(j));
        end
    end
    lls(t)  = 0.5*(likConst + sum(log(h(:,t))) + sum(data(:,t)./h(:,t)));
end
ll = sum(lls);

if isnan(ll) || isinf(ll) || ~isreal(ll)
    ll = 1e7;
end

if nargout>2
    h = h';
end