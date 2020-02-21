function [H,univariate] = dcc_fit_variance(data,p,o,q,gjrType,startingVals)
% Fits TARCH models for use in DCC and related estimators
%
% USAGE:
%   [H,UNIVARIATE] = dcc_fit_variance(DATA,P,O,Q,GJRTYPE)
%
% INPUTS:
%   DATA    - A column of mean zero data
%   P       - K by 1 vector of positive, scalar integers representing the number of symmetric innovations
%   O       - K by 1 vector of non-negative scalar integers representing the number of asymmetric innovations (0
%                    for symmetric processes)    
%   Q       - K by 1 vector of non-negative, scalar integers representing the number of lags of conditional
%                    variance (0 for ARCH) 
%   GJRTYPE - K by 1 vector of model types:
%                    1 - Model evolves in absolute values
%                    2 - Model evolves in squares [DEFAULT]
%   STARTINGVALS - [OPTIONAL] K+sum(P)+sum(O)+sum(Q) vector of starting values
%
% OUTPUTS:
%   H          - T by K matrix of conditional variances
%   UNIVARIATE - A cell array of structures used to reconstruct the variance
% 
% COMMENTS:
%
%  See also TARCH, DCC, CCC_MVGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 17/4/2012


if size(startingVals,2)>size(startingVals,1)
    startingVals = startingVals';
end

[T,k] = size(data);
H = zeros(T,k);
univariate = cell(k,1);
univariteOptions = optimset('fminunc');
univariteOptions.Display = 'none';
univariteOptions.LargeScale = 'off';
offset = 0;

for i=1:k
    if ~isempty(startingVals)
        count = 1+p(i)+o(i)+q(i);
        volStartingVals = startingVals(offset + (1:count));
        offset = offset + count;
    else
        volStartingVals = [];
    end
    [parameters, ~, ht, ~, ~, scores, diagnostics] = tarch(data(:,i),p(i),o(i),q(i), [], gjrType(i), volStartingVals, univariteOptions);
    % Store output for later use
    univariate{i}.p = p(i);
    univariate{i}.o = o(i);
    univariate{i}.q = q(i);
    univariate{i}.fdata = diagnostics.fdata;
    univariate{i}.fIdata = diagnostics.fIdata;
    univariate{i}.back_cast = diagnostics.back_cast;
    univariate{i}.m = diagnostics.m;
    univariate{i}.T = diagnostics.T;
    univariate{i}.tarch_type = gjrType(i);
    univariate{i}.parameters = parameters;
    univariate{i}.ht = ht;
    univariate{i}.A = diagnostics.A;
    univariate{i}.scores = scores;
    H(:,i) = ht;
end
