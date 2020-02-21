function ht=igarch_core(fepsilon,parameters,backCast,p,q,m,T,igarchType,constant)
% Conditional variance computation for a IGARCH(P,Q) process.
%
% USAGE:
%   [HT] = igarch_core(FEPSILON,PARAMETERS,BACKCAST,P,Q,M,T,IGARCHTYPE,CONSTANT)
%
% INPUTS:
%   FEPSILON     - Either abs(EPSILON) or EPSILON.^2, depending on IGARCHTYPE
%   PARAMETERS   - CONSTANT+P+Q by 1 vector of parameters
%   BACKCAST     - Value to be used for initializing the recursion
%   P            - Positive, scalar integer representing the number of
%                    symmetric innovations
%   Q            - Non-negative, scalar integer representing the number
%                    of lags of conditional variance (0 for ARCH)
%   M            - Number of back casts needed
%   T            - Length of FDATA, including any appended back casts
%   IGARCHTYPE   - The type of variance process, either
%                   1 - Model evolves in absolute values
%                   2 - Model evolves in squares
%   CONSTANT     - 1 if model includes a constant, 0 otherwise. 
%
% OUTPUTS:
%   HT           - Vector of conditonal varainces, T by 1
%
% COMMENTS:
%    The conditional variance, h(t), of an IGARCH(P,Q) process is modeled
%    as follows:
%
%     g(h(t)) = omega
%             + alpha(1)*f(r_{t-1}) + ... + alpha(p)*f(r_{t-p})+...
%             beta(1)*g(h(t-1)) +...+ beta(q)*g(h(t-q))
%
%     where f(x) = abs(x)  if IGARCHTYPE=1
%           g(x) = sqrt(x) if IGARCHTYPE=1
%           f(x) = x^2     if IGARCHTYPE=2
%           g(x) = x       if IGARCHTYPE=2
%     and
%     
%     sum(alpha)+sum(beta) = 1
%
%  See also IGARCH
%
%  You should use the MEX files (or compile if not using Win64 Matlab)
%  as they provide speed ups of approx 10 times relative to the m file

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009

%Initialize ht
ht=zeros(size(fepsilon));
%Set the back casts
ht(1:m)=backCast;
% Compute final beta, needed since the number of ARCH + GARCH coefficients
% is 1 less than in a usual GARCH  
archParameters = parameters(constant+1:constant+p+q-1);
finalParameter = 1 - sum(archParameters);
for i=m+1:T
    if constant
        ht(i) = parameters(1);
    end
    for j=1:p
        ht(i) = ht(i) + parameters(j+constant)*fepsilon(i-j);
    end
    for j=1:(q-1)
        ht(i) = ht(i) + parameters(j+p+constant)*ht(i-j) ;
    end
    % This line is needed since the number of ARCH + GARCH coefficients is
    % 1 less than in a usual GARCH
    ht(i) = ht(i) + finalParameter*ht(i-q) ;
end

%Return squared ht if the recursion was for a IGARCH/absolute value process
if igarchType==1
    ht=ht.^2;
end