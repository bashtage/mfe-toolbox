function ht=aparch_core(data_aug,abs_data_aug,parameters,p,o,q,m,T,back_cast,LB,UB)
% Conditional variance computation for a APARCH(P,O,Q) process.
%
% USAGE:
%   [HT] = aparch_core(DATA_AUG,ABS_DATA_AUG,PARAMETERS,P,O,Q,M,T)
%
% INPUTS:
%   DATA_AUG      - Vector of mean zero residuals augmented with zeros
%   ABS_DATA_AUG  - Absolute value of augmented data
%   PARAMETERS    - 1+P+O+Q by 1 vector of parameters
%   P             - Positive, scalar integer representing the number of
%                   symmetric innovations
%   O             - Non-negative scalar integer representing the number
%                   of asymmetric innovations (0 for symmetric processes)
%   Q             - Non-negative, scalar integer representing the number
%                   of lags of conditional variance (0 for ARCH)
%   M             - Number of back casts needed
%   T             - Length of FDATA, including any appended back casts
%
% OUTPUTS:
%   HT            - Vector of conditonal varainces, T by 1
%
% COMMENTS:
%    The conditional variance, h(t), of a APARCH(P,O,Q) process is modeled
%    as follows:
%
%
%  See also APARCH
%
%  You should use the MEX files (or compile if not using Win32 Matlab)
%  as they provide speed ups of approx 10 times relative to the m file

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

%Initialize ht
htdelta=zeros(size(data_aug));
ht=zeros(size(data_aug));
delta = parameters(1+p+o+q+1);

%Set the back casts
htdelta(1:m)=back_cast;
ht(1:m)=data_aug(m+1:T)'*data_aug(m+1:T)/(T-m);
%Recursion Loop
deltainv=2/delta;

omega = parameters(1);
alpha = parameters(2:p+1);
gamma = parameters(p+2:p+o+1);
beta = parameters(p+o+2:p+o+q+1);
for i=m+1:T
    htdelta(i) = omega;
    for j=1:p
        if o>=j
            htdelta(i) = htdelta(i) + alpha(j)*(abs_data_aug(i-j)+gamma(j)*data_aug(i-j))^delta;
        else
            htdelta(i) = htdelta(i) + alpha(j)*(abs_data_aug(i-j))^delta;
        end
    end
    for j=1:q
        htdelta(i) = htdelta(i) + beta(j)*htdelta(i-j);
    end
    % Check that HT is between its LB and UB
    htTemp = htdelta(i)^deltainv;
    if htTemp<LB
        htdelta(i)=LB^delta;
        htTemp = htdelta(i)^deltainv;
    elseif htTemp>UB
        htdelta(i)=UB^delta;
        htTemp = htdelta(i)^deltainv;
    end
    ht(i)=htTemp;
end