function rk = realized_kernel_core(returns, weights, options)
% Realized Kernel core routine that computes the value of a realized kernel
% given a set of returns and weights.
%
% USAGE:
%   [RK] = realized_kernel_core(RETURNS,WEIGHTS,OPTIONS)
%
% INPUTS:
%   RETURNS - m by 1 column vector of returns
%   WEIGHTS - H by 1 column vector of kernel weights corresponding to lags 1, 2, ..., H.  H should
%               be much smaller than m 
%   OPTIONS - A realized kernel options structure.  See help realized_options for details
%
% OUTPUTS:
%   RK        - Realized kernel value
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL
%
%  See also REALIZED_KERNEL, REALIZED_OPTIONS, REALIZED_KERNEL_WEIGHTS, REALIZED_KERNEL_BANDWIDTH
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=3
    error('Three inputs required.')
end
if size(returns,2)>size(returns,1)
    returns=returns';
end
if size(returns,2)>1
    error('RETURNS must be a m by 1 vector.')
end
if size(weights,2)>size(weights,1)
    weights=weights';
end
if size(weights,2)>1
    error('WEIGHTS must be a H by 1 vector.')
end
 
% Get the number of returns
m=size(returns,1);
if length(weights)>=m
    warning('oxfordRealized:realizedKernelLength','The length of WEIGHTS is longer than the length of RETURNS.\n  The weights are being truncated at N-1 where N is the number of returns')
    weights = weights(1:m-1);
elseif (isfield(options,'maxBandwidthPerc') && ~isempty(options.maxBandwidthPerc)) || (isfield(options,'maxBandwidth') && ~isempty(options.maxBandwidth))
    warningString = [];
    
    if ~isempty(options.maxBandwidthPerc) && (length(weights)>(options.maxBandwidthPerc * options.filteredN))
        options.maxBandwidth = round(options.maxBandwidthPerc * options.filteredN);
        warningString = ['The estimated bandwidth requires a lag length larger than ' ...
            num2str(100*options.maxBandwidthPerc) ' %% of the available data.  Bandwidth has ' ...
            'been truncated to '  num2str(options.maxBandwidth) '.'];
    elseif ~isempty(options.maxBandwidth) && length(weights)>options.maxBandwidth
        warningString = ['The estimated bandwidth requires a lag length larger than ' ...
            num2str(options.maxBandwidth) '.  Bandwidth has ' ...
            'been truncated to '  num2str(options.maxBandwidth) '.'];
    end
    
    if ~isempty(warningString)
        warning('oxfordRealized:realizedKernelLength',warningString);
        weights = weights(1:options.maxBandwidth);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Extract relevant fields form OPTIONS
isJittered = strcmpi(options.endTreatment,'jitter');
% Get the size of the kernel
H=length(weights);
if isJittered
    % Use all returns when computing autocovariances
    gammaMinus = zeros(H,1);
    gammaPlus  = zeros(H,1);
    % Loop and compute the values for gamma
    for i=1:H
        returnsMinus  = returns(1:m-i);
        returnsPlus   = returns(i+1:m);
        gammaMinus(i) = returnsMinus' * returnsPlus;
        gammaPlus(i)  = returnsMinus' * returnsPlus;
    end
    % Compute the 0-lag gamma using all returns
    gamma0 = returns' * returns;
else
    % When not jittered a different strategy is used
 
    % Construct the "central" set of returns
    returnsBase   = returns(H+1:m-H);
    % Initialized placeholders for the terms of gamma
    gammaMinus = zeros(H,1);
    gammaPlus  = zeros(H,1);
    % Loop and compute the values for gamma
    for i=1:H
        returnsMinus  = returns(H+1-i:m-H-i);
        returnsPlus   = returns(H+1+i:m-H+i);
        gammaMinus(i) = returnsMinus' * returnsBase;
        gammaPlus(i)  = returnsBase'  * returnsPlus;
    end
    % Compute the 0-lag gamma
    gamma0 = returnsBase' * returnsBase;
end
 
% Construct the kernel
if ~isempty(weights)
    rk = gamma0 + weights'*(gammaMinus+gammaPlus);
else
    warning('oxfordRealized:realizedKernelLength','The number of lags used was 0.');
    rk = gamma0;
end