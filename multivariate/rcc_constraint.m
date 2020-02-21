function [c,ceq] = rcc_constraint(parameters,data,m,n,R,backCast,stage,type,composite,isJoint,isInference,rScale,univariate)  %#ok<INUSD,INUSL>
% Non-linear constraint for estimation of RARCH(p,q) multivariate volatility models
%
% USAGE:
%  [C,CEQ] = rarch_constraint(PARAMETERS,DATA,P,Q,C,BACKCAST,TYPE,ISJOINT)
%
% INPUTS:
%   See rarch_likelihood
%
% OUTPUTS:
%   C   - Vector of inequality constraints
%   CEQ - Empty
%
% COMMENTS:
%
% See also RARCH, RARCH_LIKELIHOOD

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 3/27/2012

ceq = [];
k = size(data,1);
if stage==3
    tempParameters = parameters;
else
    tempParameters = parameters(k*(k-1)/2:length(parameters));
end
[~,A,B] = rarch_parameter_transform(tempParameters,m,n,k,[],type,false,false);
constraint = diag(sum(A.^2,3)+sum(B.^2,3) - .99998);
switch type
    case 1
        c = constraint(1);
    case 2
        theta = parameters(length(parameters)).^2;
        c = diag(sum(A.^2,3)) - theta;
    case 3
        c = constraint;
end
