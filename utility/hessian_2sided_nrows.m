function H = hessian_2sided_nrows(f,x,k,varargin)
% Computes the last K rows of a 2-sided finite difference Hessian
%
% USAGE:
%   H=hessian_2sided_nrows(FUNC,X,K,VARARGIN)
%
% INPUTS:
%   FUNC         - Function name, fval = func(x,varargin)
%   X            - Vector of parameters (N x 1)
%   K            - Number of rows to compute
%   VARARGIN     - Optional arguments passed to the function
%
% OUTPUTS:
%   H            - Finite difference, 2-sided hessian: N x K
%
% COMMENTS:
%   Modification of hessian_2sided
%   See also HESSIAN_2SIDED

% Code originally from COMPECON toolbox [www4.ncsu.edu/~pfackler]
% documentation modified to fit the format of the Econometrics Toolbox
% by James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com
%
% Further modified (to do 2-sided numerical derivs, rather than 1) by:
% Modifications Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

if size(x,2)>size(x,1)
    x=x';
end
if size(x,2)~=1
    error('X must be a column vector.')
end
n = size(x,1);

try
    feval(f,x,varargin{:});
catch FE
    errtxt=['There was an error evaluating the function.  Please check the arguements.  The specific error was:' FE.message];
    error(errtxt);
end


fx = feval(f,x,varargin{:});

% Compute the stepsize (h)
h = eps.^(1/3)*max(abs(x),1e-2);
xh = x+h;
h = xh-x;

ee = sparse(1:n,1:n,h,n,n);

% Compute forward and backward steps
gp = zeros(n,1);
gm = zeros(n,1);
for i=1:n
    gp(i) = feval(f,x+ee(:,i),varargin{:});
    gm(i) = feval(f,x-ee(:,i),varargin{:});
end

hh=h*h';
Hm=NaN*ones(n);
Hp=NaN*ones(n);
% Compute "double" forward and backward steps
for i=n-k+1:n
    for j=1:n
        Hp(i,j) = feval(f,x+ee(:,i)+ee(:,j),varargin{:});
        Hp(j,i)=Hp(i,j);
        Hm(i,j) = feval(f,x-ee(:,i)-ee(:,j),varargin{:});
        Hm(j,i)=Hm(i,j);
    end
end
%Compute the hessian
H = zeros(n);
for i=(n-k+1):n
    for j=1:n
        H(i,j) = (Hp(i,j)-gp(i)-gp(j)+fx+fx-gm(i)-gm(j)+Hm(i,j))/hh(i,j)/2;
        H(j,i) = H(i,j);
    end
end

%Only the final K rows
H=H((n-k+1):n,:);
