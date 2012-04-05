function ht=tarch_core(fdata,fIdata,parameters,back_cast,p,o,q,m,T,tarch_type)
% Conditional variance computation for a TARCH(P,O,Q) process.
%
% USAGE:
%   [HT] = tarch_core(FDATA,FIDATA,PARAMETERS,BACK_CAST,P,O,Q,M,T,TARCH_TYPE)
%
% INPUTS:
%   FDATA         - A column of mean zero data transformed according to F (see COMMENTS)
%   FIDATA        - FDATA * DATA<0
%   PARAMETERS    - 1+P+O+Q by 1 vector of parameters
%   BACK_CAST     - Value to be used for initializing the recursion
%   P             - Positive, scalar integer representing the number of
%                   symmetric innovations
%   O             - Non-negative scalar integer representing the number
%                   of asymmetric innovations (0 for symmetric processes)
%   Q             - Non-negative, scalar integer representing the number
%                   of lags of conditional variance (0 for ARCH)
%   M             - Number of back casts needed
%   T             - Length of FDATA, including any appended back casts
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
%  See also TARCH
%
%  You should use the MEX files (or compile if not using Win32 Matlab)
%  as they provide speed ups of approx 10 times relative to the m file

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

%Initialize ht
ht=zeros(size(fdata));
%Set the back casts
ht(1:m)=back_cast;

%Recursion Loop
for i=m+1:T
    ht(i) = parameters(1);
    for j=1:p
        ht(i) = ht(i) + parameters(j+1)*fdata(i-j);
    end
    for j=1:o
        ht(i) = ht(i) + parameters(j+p+1)*fIdata(i-j) ;
    end
    for j=1:q
        ht(i) = ht(i) + parameters(j+p+o+1)*ht(i-j) ;
    end
end


%Return squared ht if the recursion was for a TARCH/absolute value process
if tarch_type==1
    ht=ht.^2;
end
