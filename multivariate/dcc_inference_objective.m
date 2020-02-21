function [obj,objs] = dcc_inference_objective(parameters,data,dataAsym,m,l,n,univariate) %#ok<INUSL>
% Objective function used by DCC, CCC and related multivariate volatility models to perform
% inference.  It is not used (and cannot be used) to estimate parameters.
%
% USAGE:
%  [OBJ,OBJS] = dcc_inference_objective(PARAMETERS,DATA,DATAASYM,M,L,N,UNIVARIAT)
%
% INPUTS:
%   PARAMETERS   - Vector of ADCC parameters including possibly volatility and intercepts
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   DATAASYM     - [OPTIONAL] K by K by T array of asymmetric covariance estimators only needed if
%                    DATA is 3-dimensional and O>0 or L>0
%   M            - Order of symmetric innovations in DCC model
%   L            - Order of asymmetric innovations in ADCC model
%   N            - Order of lagged correlation in DCC model
%   UNIVARIATE   - Cell array of structures containing information needed to compute volatilities.
%
% OUTPUTS:
%   OBJ          - "As if" objective for estimating correlation intercept parameters
%   OBJS         - T by 1 vector of "as if" objectives for estimating correlation intercept parameters
%
% COMMENTS:
%
% See also DCC, CCC_MVGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/13/2012


[k,~,T] = size(data);
% Parse parameters
count = 0;
for i=1:k
    u = univariate{i};
    count = count + u.p+u.o+u.q+1;
end
garchParameters = parameters(1:count);
offset = count;
% R is next
count = k*(k-1)/2;
R = corr_ivech(parameters(offset + (1:count)));
offset = offset + count;
if l>0
    count = k*(k+1)/2;
    N = ivech(parameters(offset+(1:count)));
end

H = dcc_reconstruct_variance(garchParameters,univariate);
stdData = zeros(k,k,T);
stdDataAsym = zeros(k,k,T);
for t=1:T
    h = sqrt(H(t,:));
    stdData(:,:,t) = data(:,:,t)./(h'*h);
    stdDataAsym(:,:,t) = dataAsym(:,:,t)./(h'*h);
end

scales = diag(mean(stdData,3));
objs = zeros(T,1);
for j=1:k-1 % Cols
    for i=j+1:k % Rows
        scale = sqrt(scales(i)*scales(j));
        errors = squeeze(stdData(i,j,:))/scale - R(i,j);
        objs = objs + 0.5*(errors.^2);
    end
end

if l>0
    for j=1:k
        for i= j:k
            errors = squeeze(stdDataAsym(i,j,:)) - N(i,j);
            objs = objs + 0.5*(errors.^2);
        end
    end
end

obj = sum(objs);

