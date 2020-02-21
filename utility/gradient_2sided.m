function [G,Gt] = gradient_2sided(f,x,varargin)
% Computes 2-sided finite difference gradient vector
%
% USAGE:  
%   [G,GT]=gradient_2sided(FUNC,X,VARARGIN)
%
% INPUTS:
%   FUNC         - function name, fval = func(x,varargin)
%   X            - vector of parameters (M x 1)
%   VARARGIN     - optional arguments passed to the function
%
% OUTPUTS:
%   G            - finite differnce, 2-sided gradient
%   GT           - T by M matrix of individual scores.  Used in estimating
%                   covariances
%
% COMMENTS:
%   Modification of hessian() from the jpl toolbox


% Code originally from COMPECON toolbox [www4.ncsu.edu/~pfackler]
% documentation modified to fit the format of the Ecoometrics Toolbox
% by James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com
%
% Further modified (to do 2-sided numerical derivs, rather than 1)
%
% Modifications Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 2/1/2006

if size(x,2)>size(x,1)
    x=x';
end

M = size(x,1);
 
h = eps.^(1/3)*max(abs(x),1e-2);
xh = x+h;
h = xh-x;    
% Compute forward step 
ee = diag(h);%sparse(1:M,1:M,h,M,M)
gf = zeros(M,1);
gb = zeros(M,1);

if nargout==2
    [temp1,temp2] = feval(f,x+ee(:,1),varargin{:});
    T=length(temp2);
    Gf=zeros(T,M);
    Gb=zeros(T,M);
end

for i=1:M
    if nargout<2
  gf(i) = feval(f,x+ee(:,i),varargin{:});
    else
        [gf(i),Gf(:,i)] = feval(f,x+ee(:,i),varargin{:});
    end
end

% Compute backward step 
for i=1:M
    if nargout<2
  gb(i) = feval(f,x-ee(:,i),varargin{:});
    else
        [gb(i),Gb(:,i)] = feval(f,x-ee(:,i),varargin{:});
    end
end

G = (gf-gb)./(2*h);
if nargout==2
    T=size(Gb,1);
    Gt = (Gf-Gb)./repmat(2*h',T,1);
end
