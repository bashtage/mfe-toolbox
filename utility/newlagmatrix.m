function [y,x]=newlagmatrix(x,nlags,c)
% Construction of a matrix of lags (X) and a vector (Y) for use in an
% autoregression
%
% USAGE:
%     [Y,X]=newlagmatrix(Y,NLAGS,C)
%
% INPUTS:
%     Y     - The dependant variable(Tx1)
%     NLAGS - The number of lags(scalar)
%     C     - [OPTIONAL] 1 if you want to include a constant; 0 is default
%
% OUTPUTS:
%     Y     - (n-p) by 1 vector of contemporaneous values
%     X     - (n-p) by p (or p+1 if c=1) matrix of lags and possibly constant.  The
%               matrix is of the form [constant(if included) t-1 t-2 ... t-p]
%
% COMMENTS:
%     Original name, 'lagmatrix' conflicts with a Matlab file, so newlagmatrix
%     replaces original 

% Author: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 5    Date: 12/1/2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 2
        c=0;
    case 3
        %Nothing
    otherwise
        error('newlagmatrix requires either 2 or 3 arguments.')
end

[T,K]=size(x);
if min(T,K)~=1
    error('X must be a vector')
end

if T<K
    x=x'; %Transpose if a column vector is passed
end

if  length(c)~=1 || any(c~=0 && any(c~=1)) 
    error('C must be 1 or 0')
end

if length(nlags)~=1 || floor(nlags)~=nlags || any(nlags<0)
    error('NLAGS must be a positive integer')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nlags>0
    nlags=nlags+1;
    newX=[x;zeros(nlags,1)];
    lagmatrix=repmat(newX,nlags,1);
    lagmatrix=reshape(lagmatrix(1:size(lagmatrix,1)-nlags),T+nlags-1,nlags);
    lagmatrix=lagmatrix(nlags:T,:);
    y=lagmatrix(:,1);
    x=lagmatrix(:,2:nlags);
    if c==1
        x=[ones(size(x,1),1) x];
    end
else
    if c==1
        y=x;
        x=ones(T,1);
    else
        y=x;
        x=[];
    end
end
