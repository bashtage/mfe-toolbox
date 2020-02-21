function [arroots, absarroots] = armaroots(parameters, constant, p, q, X)
% Computes the roots of the characteristic equation of an ARMAX(P,Q) as parameterized by ARMAXFILTER
%
% USAGE:
%   [ARROOTS,ABSARROOTS] = armaroots(PARAMETERS,CONSTANT,P,Q)
%   [ARROOTS,ABSARROOTS] = armaroots(PARAMETERS,CONSTANT,P,Q,X)
% 
% INPUTS:
%   PARAMETERS - A CONSTANT+length(P)+length(Q)+size(X,2) by 1 vector of parameters, usually an
%                  output from ARMAXFILTER 
%   CONSTANT   - Scalar variable: 1 to include a constant, 0 to exclude
%   P          - Non-negative integer vector representing the AR orders to include in the model.
%   Q          - Non-negative integer vector representing the MA orders to include in the model.
%   X          - [OPTIONAL] A T by K matrix of exogenous variables. 
%
% OUTPUTS:
%   ARROOTS     - A max(P) by 1 vector containing the roots of the characteristic equation
%                   corresponding to the ARMA model input 
%   ABSARROOTS  - Absolute value or complex modulus of the autoregressive roots
% 
% COMMENTS:
%
% EXAMPLES:
%   Compute the AR roots of an ARMA(2,2)
%       phi = [1.3 -.35]; theta = [.4 .3]; parameters=[1 phi theta]'; 
%       [arroots, absarroots] = armaroots(parameters, 1, [1 2], [1 2])
%   Compute the AR roots of an irregular AR(3)
%       phi = [1.3 -.35]; parameters = [1 phi]';                    
%       [arroots, absarroots] = armaroots(parameters, 1, [1 3],[])  
%
% See also ARMAXFILTER, ROOTS
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 1/1/2007
 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>5 || nargin<4
    error('4 or 5 inputs required.')
end
if nargin==4
    X=[];
end
%%%%%%%%%%%%%%%
% P
%%%%%%%%%%%%%%%
if size(p,2)>size(p,1)
    p=p';
end
if isempty(p)
    p=0;
end
if min(size(p))~=1
    error('P must be a column vector of included lags')
end
if  any(p<0) || any(floor(p)~=p)
    error('P must contain non-negative integers only')
end
if size(p,1)==1 && p==0
    p=[];
end
if length(unique(q))~=length(q)
    error('P must contain at most one of each lag')
end
lp = length(p);
%%%%%%%%%%%%%%%
% Q
%%%%%%%%%%%%%%%
if size(q,2)>size(q,1)
    q=q';
end
if isempty(q)
    q=0;
end
if min(size(q))~=1
    error('Q must be a column vector of included lags')
end
if  any(q<0) || any(floor(q)~=q)
    error('Q must contain non-negative integers only')
end
if size(q,1)==1 && q==0
    q=[];
end
if length(unique(q))~=length(q)
    error('Q must contain at most one of each lag')
end
lq = length(q);
%%%%%%%%%%%%%%%
% Constant
%%%%%%%%%%%%%%%
if ~ismember(constant,[0 1])
    error('CONSTANT must be 0 or 1')
end
%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%
if size(parameters,2)>size(parameters,1)
    parameters=parameters';
end
if size(parameters,2)~=1
    error('PARAMETERS must be a column vector');
end
if length(parameters)~=lp+lq+size(X,2)+constant
    error('PARAMETERS must have length compatible with P, Q, CONSTANT and X')
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
if lp==0
    arroots=[];
    absarroots=[];
else
    if constant
        arparameters=parameters(2:lp+1);
    else
        arparameters=parameters(1:lp);
    end
    formatted_arparameters = zeros(1,lp);
    formatted_arparameters(p) = arparameters;
    arroots=roots([1 -formatted_arparameters]);
    absarroots=abs(arroots);
end
