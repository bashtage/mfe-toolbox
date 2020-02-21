function [p,q,model_type,error_type,startingvals,options]=agarch_parameter_check(epsilon, p, q, model_type, error_type, startingvals, options)
% AGARCH(P,Q) and NAGARCH(P,Q) input validation.  Ensures that
% the input parameters are conformable to what is expected.
%
% USAGE:
%   [P,Q,ERROR_TYPE,MODEL_TYPE,STARTINGVALS,OPTIONS] = ...
%        agarch_parameter_check(EPSILON,P,Q,MODEL_TYPE,ERROR_TYPE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   See AGARCH.
%
% OUTPUTS:
%   See AGARCH.
%
% COMMENTS:
%   See also AGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009



%%%%%%%%%%%%%%%
% epsilon
%%%%%%%%%%%%%%%
if size(epsilon,2) > 1 || length(epsilon)==1
    error('epsilon series must be a column vector.')
elseif isempty(epsilon)
    error('epsilon is empty.')
end

%%%%%%%%%%%%%%%
% q
%%%%%%%%%%%%%%%
if (length(q) > 1) || any(q < 0) || isempty(q)
    error('q must be a non-negative scalar.')
end

%%%%%%%%%%%%%%%
% p
%%%%%%%%%%%%%%%
if (length(p) > 1) || any(p <  1) || isempty(p)
    error('p must be positive scalar.')
end

%%%%%%%%%%%%%%%
% model_type
%%%%%%%%%%%%%%%
if nargin<4
    model_type = 'AGARCH';
end

if isempty(model_type)
    model_type = 'AGARCH';
end

model_type = upper(model_type);
switch model_type
    case 'AGARCH'
        model_type = 1;
    case 'NAGARCH'
        model_type = 2;
    otherwise
        error('MODEL_TYPE must be either ''AGARCH'' or ''NAGARCH''')
end


%%%%%%%%%%%%%%%
% error_type
%%%%%%%%%%%%%%%
if nargin<5
    error_type='NORMAL';
end

if isempty(error_type)
    error_type='NORMAL';
end

switch error_type
    case 'NORMAL'
        error_type = 1;
    case 'STUDENTST'
        error_type = 2;
    case 'GED'
        error_type = 3;
    case 'SKEWT'
        error_type = 4;
    otherwise
        error('error_type must be a string and one of: ''NORMAL'', ''STUDENTST'', ''GED'' or ''SKEWT''.');
end


%%%%%%%%%%%%%%%%%
% starting values
%%%%%%%%%%%%%%%%%
if nargin>5
    if ~isempty(startingvals)
        %Validate starting vals, different if normal than if T or GED
        if error_type==1
            if length(startingvals)~=(p+q+2) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+q+2 elements');
            end
        elseif error_type==2 %T
            if length(startingvals)~=(p+q+3) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+q+3 elements');
            end
            if startingvals(p+q+3)<2.1
                error('Nu must be greater than 2.1 when using Students-T errors');
            end
        elseif error_type==3 %GED
            if length(startingvals)~=(p+q+3) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+q+3 elements');
            end
            if startingvals(p+q+3)<1.05
                error('Nu must be greater than 1 when using GED errors');
            end
        elseif error_type==4 %SKEWT
            if length(startingvals)~=(p+q+4) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+q+4 elements');
            end
            if startingvals(p+q+3)<2.1
                error('Nu must be greater than 2.1 when using Skew T errors');
            end
            if startingvals(p+q+4)<-.9 || startingvals(p+q+4)>.9
                error('Lambda must be between -.9 and .9 when using Skew T errors');
            end
        end
        if any(startingvals(1:p+1)<=0) || any(startingvals(p+3:q+p+2)<=0)
            error('All startingvals for omega, alpha and beta must be strictly greater than zero');
        end
        alpha = startingvals(2:p+1);
        beta  = startingvals(p+3:p+q+2);
        if (sum(alpha)+sum(beta))>=1
            error('The sum of the alphas and betas must be less than 1');
        end
        gamma = startingvals(p+2);
        if model_type == 1
            if gamma<=quantile(epsilon,.01) || gamma>=quantile(epsilon,.99)
                error('The starting value for gamma in AGARCH must be between q(.01,EPSILON) and q(.99,EPSILON)');
            end
        else
            if gamma<=-4 || gamma>=4
                error('The starting value for gamma in NAGARCH must be between -4 and 4');
            end
        end
    end
else
    startingvals=[];
end

%%%%%%%%%%%%%%%
% options
%%%%%%%%%%%%%%%
if nargin>6 && ~isempty(options)
    try
        optimset(options);
    catch OE
        error(['OPTIONS is not a valid minimization option structure.  The error message was:' OE.message]);
    end
else
    %Setup the options in case of none provided
    options  =  optimset('fminunc');
    options.TolFun=1e-005;
    options.TolX=1e-005;
    options.Display='iter';
    options.Diagnostics='on';
    options.LargeScale='off';
    options.MaxFunEvals=200*(2+p+q);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(q)
    q=0;
end
