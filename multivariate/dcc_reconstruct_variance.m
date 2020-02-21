function H = dcc_reconstruct_variance(garchParameters,univariate)
% Computes the variance according to TARCH models for use in DCC and related estimators
%
% USAGE:
%   [H] = dcc_reconstruct_variance(GARCHPARAMETERS,UNIVARIATE)
%
% INPUTS:
%   PARAMETERS - A vector of parameters with K+sum(P)+sum(O)+sum(Q) elemments in the order 
%                  [vol1 vol2 ... volK]
%   UNIVARIATE - A cell array of structures used to reconstruct the variance
%
% OUTPUTS:
%   H          - T by K matrix of conditional variances
% 
% COMMENTS:
%
%  See also TARCH, DCC, CCC_MVGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 17/4/2012

k = length(univariate);
T = univariate{1}.T-univariate{1}.m;
H = zeros(T,k);

offset = 0;
for i=1:k
    u = univariate{i};
    count = u.p+u.o+u.q+1;
    volParameters = garchParameters(offset + (1:count));
    offset = offset+count;
    ht = tarch_core(u.fdata,u.fIdata,volParameters,u.back_cast,u.p,u.o,u.q,u.m,u.T,u.tarch_type);
    H(:,i) = ht(u.m+1:u.T);
end