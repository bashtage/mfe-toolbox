function ht=tarch_core_simple(data,parameters,backCast,backCastAsym,p,o,q,tarch_type)
% Simple (but clow) implementation of conditional variance computation for a TARCH(P,O,Q) process.
%
% USAGE:
%   [HT] = tarch_core(DATA,PARAMETERS,BACKCAST,BACKCASTASYM,P,O,Q,TARCH_TYPE)
%
% INPUTS:
%   DATA          - A column of mean zero data
%   PARAMETERS    - 1+P+O+Q by 1 vector of parameters
%   BACKCAST      - Value to be used for initializing the recursion
%   BACKCASTASYM  - Value to be used for initializing the asymmetric term in the recursion
%   P             - Positive, scalar integer representing the number of symmetric innovations
%   O             - Non-negative scalar integer representing the number of asymmetric innovations
%   Q             - Non-negative, scalar integer representing the number of lags of conditional variance
%   TARCH_TYPE    - The type of variance process, either
%                     1 - Model evolves in absolute values
%                     2 - Model evolves in squares
%
% OUTPUTS:
%   HT            - Vector of conditonal varainces, T by 1
%
% COMMENTS:
%    The conditional variance, h(t), of a TARCH(P,O,Q) process is modeled
%    as follows:
%
%     g(h(t)) = omega
%             + alpha(1)*f(r_{t-1}) + ... + alpha(p)*f(r_{t-p})+...
%             + gamma(1)*I(t-1)*f(r_{t-1}) +...+ gamma(o)*I(t-o)*f(r_{t-o})+...
%             beta(1)*g(h(t-1)) +...+ beta(q)*g(h(t-q))
%
%     where f(x) = abs(x)  if tarch_type=1
%           g(x) = sqrt(x) if tarch_type=1
%           f(x) = x^2     if tarch_type=2
%           g(x) = x       if tarch_type=2
%
%  See also TARCH, TARCH_CORE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 4    Date: 4/12/2012

%Initialize ht
T = length(data);
ht = zeros(T,1);
if tarch_type==1
    fdata = abs(data);
else
    fdata = data.^2;
end
fIdata = fdata.*(data<0);

%Recursion Loop
for i=1:T
    ht(i) = parameters(1);
    for j=1:p
        if (i-j)>0
            ht(i) = ht(i) + parameters(j+1)*fdata(i-j);
        else
            ht(i) = ht(i) + parameters(j+1)*backCast;
        end
    end
    for j=1:o
        if (i-j)>0
            ht(i) = ht(i) + parameters(j+p+1)*fIdata(i-j);
        else
            ht(i) = ht(i) + parameters(j+p+1)*backCastAsym;
        end
    end
    for j=1:q
        if (i-j)>0
            ht(i) = ht(i) + parameters(j+p+o+1)*ht(i-j) ;
        else
            ht(i) = ht(i) + parameters(j+p+o+1)*backCast ;
        end
    end
end

% Return squared ht if the recursion was for a TARCH/absolute value process
if tarch_type==1
    ht=ht.^2;
end