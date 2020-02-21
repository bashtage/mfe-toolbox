function [sigma,correl] = cov2corr(cov)
% Convert covariance matrices to correlation matrices.  
%
% USAGE:
%  [SIGMA,CORREL] = cov2corr(COV)
%
% INPUTS:
%   COV   - A K by K covariance matrix -OR- 
%           A K by K by T array of covariance matrices
%
% OUTPUTS:
%   SIGMA  - A K by 1 vector of standard deviations if COV is K by K -OR-
%            A T by K matrix of standard deviations
%   CORREL - A K by K matrix of correlations -OR-
%            A K by K by T matrix of correlations.
%
% EXAMPLES:
%   % DCC(1,1)
%   [~,~,Ht] = dcc(data,[],1,0,1)
%   [S,Rt] = cov2corr(Rt);
%
% See also DCC

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 10/23/2012

if ndims(cov)==2
    sigma = sqrt(diag(cov));
    correl = cov ./(sigma*sigma');
elseif ndims(cov)==3
    T = size(cov,3);
    K = size(cov,2);
    sigma = zeros(T,K);
    correl = zeros(K,K,T);
    for t=1:T
        sigma(t,:) = sqrt(diag(cov(:,:,t)))';
        correl(:,:,t) = cov(:,:,t)./(sigma(t,:)'*sigma(t,:));
    end
end