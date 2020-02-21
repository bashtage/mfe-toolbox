function ht=agarch_core(data,parameters,back_cast,p,q,m,T,model_type)
% Conditional variance computation for a AGARCH(P,Q) process.
%
% USAGE:
%   [HT] = agarch_core(DATA,PARAMETERS,BACK_CAST,P,Q,M,T,MODEL_TYPE)
%
% INPUTS:
%   DATA          - A column of mean zero data audmented with M backcasts
%   PARAMETERS    - 2+P+Q by 1 vector of parameters
%   BACK_CAST     - Value to be used for initializing the recursion
%   P             - Positive, scalar integer representing the number of
%                   symmetric innovations
%   Q             - Non-negative, scalar integer representing the number
%                   of lags of conditional variance (0 for ARCH)
%   M             - Number of backcasts needed
%   T             - Length of DATA, including any appended back casts
%   MODEL_TYPE    - The type of variance process, either
%                     1 - AGARCH
%                     2 - NAGARCH
%
% OUTPUTS:
%   HT            - Vector of conditonal varainces, T by 1
%
% COMMENTS:
%    The conditional variance, h(t), of a AGARCH(P,Q) process is given by:
%
%     h(t)  = omega
%             + alpha(1)*(r_{t-1}-gamma)^2 + ... + alpha(p)*(r_{t-p}-gamma)^2
%             + beta(1)*h(t-1) +...+ beta(q)*h(t-q)
%
%    The conditional variance, h(t), of a NAGARCH(P,Q) process is given by:
%
%     h(t)  = omega
%             + alpha(1)*(r_{t-1}-gamma*sqrt(h(t-1)))^2 + ... + alpha(p)*(r_{t-p}-gamma*sqrt(h(t-p)))^2 
%             + beta(1)*h(t-1) +...+ beta(q)*h(t-q)
%
%  See also AGARCH
%
%  You should use the MEX files (or compile if not using Win32 Matlab)
%  as they provide speed ups of approx 10 times relative to the m file

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009


%Initialize ht
ht=zeros(size(data));
%Set the back casts
ht(1:m)=back_cast;
gamma = parameters(2+p);
shock = zeros(m,1);

if model_type == 1
    for i = 1:m
        shock(i) = back_cast + gamma^2;
    end
    for i=m+1:T
        ht(i) = parameters(1);
        for j=1:p
            ht(i) = ht(i) + parameters(j+1)*shock(i-j);
        end
        for j=1:q
            ht(i) = ht(i)+ parameters(j+2+p)*ht(i-j);
        end
        shock(i) = (data(i)-gamma)^2;
    end
else
    for i = 1:m
        shock(i) = back_cast * (1+ gamma^2);
    end
    for i=m+1:T
        ht(i) = parameters(1);
        for j=1:p
            ht(i) = ht(i) + parameters(j+1)*shock(i-j);
        end
        for j=1:q
            ht(i) = ht(i)+ parameters(j+2+p)*ht(i-j);
        end
        shock(i) = (data(i)-gamma*sqrt(ht(i)))^2;
    end
end
