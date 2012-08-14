function [bsdata, indices]=stationary_bootstrap(data,B,w)
% Implements the stationary bootstrap for bootstrapping stationary, dependent series
% 
% USAGE:
%   [BSDATA, INDICES] = stationary_bootstrap(DATA,B,W)
% 
% INPUTS:
%   DATA   - T by 1 vector of data to be bootstrapped
%   B      - Number of bootstraps
%   W      - Average block length. P, the probability of starting a new block is defined P=1/W
% 
% OUTPUTS:
%   BSDATA  - T by B matrix of bootstrapped data
%   INDICES - T by B matrix of locations of the original BSDATA=DATA(indexes);
% 
% COMMENTS:
%   To generate bootstrap sequences for other uses, such as bootstrapping vector processes,
%   set DATA to (1:N)'.  
%
% See also block_bootstrap

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 2    Date: 12/31/2001



%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=3
    error('3 inputs required')
end
% Get the length of the data
[t,k]=size(data);
if k>1
    error('DATA must be a column vector')
end
if t<2 
    error('DATA must have at least 2 observations.')
end
if ~isscalar(w) || w<=0
    error('W must be a positive scalar.')
end
if ~isscalar(B) || B<1 || floor(B)~=B
    error('B must be a positive scalar integer')
end
%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the probability of a new block
p=1/w;
% Set up the bsdata and indices
indices=zeros(t,B);
% Initial positions
indices(1,:)=ceil(t*rand(1,B));
% Set up the random numbers
select=rand(t,B)<p;
indices(select)=ceil(rand(1,sum(sum(select)))*t);
for i=2:t
    % Determine whether we stay (rand>p) or move to a new starting value
    % (rand<p)
    indices(i,~select(i,:))=indices(i-1,~select(i,:))+1;
end
indices(indices>t) = indices(indices>t)-t;
% The indices make finding the bsdata simple
bsdata=data(indices);