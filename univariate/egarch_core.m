function ht=egarch_core(data,parameters,back_cast,upper,p,o,q,m,T)
% Conditional variance computation for a EGARCH(P,O,Q) process.
%
% USAGE:
%   [HT] = egarch_core(DATA,PARAMETERS,BACK_CAST,UPPER,P,O,Q,M,T)
%
% INPUTS:
%   DATA          - A column of mean zero data 
%   PARAMETERS    - 1+P+O+Q by 1 vector of parameters
%   BACK_CAST     - Value to be used for initializing the recursion
%   UPPER         - Upper bound for ht, larger values will be set to this
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
%   The conditional variance, h(t), of an EGARCH(P,O,Q) process is modeled
%   as follows:
%
%   ln(h(t)) = omega
%            + alpha(1)*(abs(e_{t-1})-C) + ... + alpha(p)*(abs(e_{t-p})-C)+...
%            + gamma(1)*e_{t-1} +...+ e_{t-o} +...
%              beta(1)*ln(h(t-1)) +...+ beta(q)*ln(h(t-q))
%
%  See also EGARCH
%
%  You should use the MEX files (or compile if not using Win32 Matlab) 
%  as they provide speed ups of approx 100 times relative to the m file

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005



ht=zeros(size(data));
eht=zeros(size(data));
ht(1:m)=back_cast;
eht(1:m)=exp(back_cast);
data2=data;
absdata2=data;
subconst=0.797884560802865;
for i=1:m;
    vol=sqrt(eht(i));
    data2(i)=data(i)/vol;
    absdata2(i)=abs(data2(i))-subconst;
end

eLB=eht(1)/10000;
hLB=ht(1)-log(10000);

for i=m+1:T
    ht(i) = parameters(1);
    for j=1:p
        ht(i) = ht(i) + parameters(j+1)*absdata2(i-j);
    end
    for j=1:o
        ht(i) = ht(i) + parameters(j+p+1)*data2(i-j) ;
    end
    for j=1:q
        ht(i) = ht(i) + parameters(j+p+o+1)*ht(i-j) ;
    end
    %ht(i) = parameters'*[1 ; absdata2(i-1:-1:i-p) ; data2(i-1:-1:i-o) ; ht(i-1:-1:i-q)];
    eht(i)=exp(ht(i));
    if eht(i)< eLB;
        eht(i) = eLB;
        ht(i) = hLB;
    end
    if eht(i)>upper;
        if ~isinf(eht(i))
            eht(i)= upper+ht(i);
        else
            eht(i)=upper+800;
        end
        %Keep ht(i) from diverging
        ht(i) = log(eht(i)) ;
    end
    vol=sqrt(eht(i));
    data2(i)=data(i)/vol;
    absdata2(i)=abs(data2(i))-subconst;
end
ht=eht;