function [p,q,errorType,igarchType,constant,startingvals,options]=igarch_parameter_check(epsilon, p, q, errorType, igarchType, constant, startingvals, options)
% IGARCH(P,Q) input validation.  Ensures that the input parameters are
% conformable to what is expected.
%
% USAGE:
%   [P,Q,ERRORTYPE,IGARCHTYPE,STARTINGVALS,OPTIONS] = ...
%        igarch_parameter_check(EPSILON,P,Q,ERRORTYPE,IGARCHTYPE,CONSTANT,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   See IGARCH.
%
% OUTPUTS:
%   See IGARCH.
%
% COMMENTS:
%   See also IGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009


%%%%%%%%%%%%%%%
% epsilon
%%%%%%%%%%%%%%%
if size(epsilon,2) > 1 || length(epsilon)==1
    error('EPSILON series must be a column vector.')
elseif isempty(epsilon)
    error('EPSILON is empty.')
end

%%%%%%%%%%%%%%%
% q
%%%%%%%%%%%%%%%
if (length(q) > 1) || any(q < 1) || isempty(q)
    error('Q must be a positive scalar.')
end


%%%%%%%%%%%%%%%
% p
%%%%%%%%%%%%%%%
if (length(p) > 1) || any(p <  1) || isempty(p)
    error('P must be positive scalar.')
end

%%%%%%%%%%%%%%%
% ERRORTYPE
%%%%%%%%%%%%%%%
if nargin<4
    errorType='NORMAL';
end

if isempty(errorType)
    errorType='NORMAL';
end

if strcmp(errorType,'NORMAL')
    errorType = 1;
elseif strcmp(errorType,'STUDENTST')
    errorType = 2;
elseif strcmp(errorType,'GED')
    errorType = 3;
elseif strcmp(errorType,'SKEWT')
    errorType = 4;
else
    error('errorType must be a string and one of: ''NORMAL'', ''STUDENTST'', ''GED'' or ''SKEWT''.');
end

%%%%%%%%%%%%%%%
% IGARCHTYPE
%%%%%%%%%%%%%%%
if nargin<5 || isempty(igarchType)
    igarchType = 2;
end

if length(igarchType)>1
    error('IGARCHTYPE must be a scalar');
end

if ~(igarchType==2 || igarchType==1)
    error('IGARCHTYPE must be either 1 or 2')
end


%%%%%%%%%%%%%%%%%
% Constant
%%%%%%%%%%%%%%%%%
if nargin<6 || isempty(constant)
    constant = 1;
end
constant = double(constant);
if ~isscalar(constant) || ~ismember(constant,[0 1]) 
    error('CONSTANT must be either 0 or 1.')
end

%%%%%%%%%%%%%%%%%
% starting values
%%%%%%%%%%%%%%%%%
if nargin>=7 && ~isempty(startingvals)
    %Validate starting vals, different if normal than if T or GED
    if errorType==1
        if length(startingvals)~=(p+q-1+constant) || size(startingvals,2)~=1
            error('startingvals must be a column vector with p+q+CONSTANT-1 elements');
        end
    elseif errorType==2 %T
        if length(startingvals)~=(p+q+constant) || size(startingvals,2)~=1
            error('startingvals must be a column vector with p+q+CONSTANT elements');
        end
        if startingvals(p+q+constant)<2.1
            error('Nu must be greater than 2.1 when using Students-T errors');
        end
    elseif errorType==3 %GED
        if length(startingvals)~=(p+q+constant) || size(startingvals,2)~=1
            error('startingvals must be a column vector with p+q+CONSTANT elements');
        end
        if startingvals(p+q+constant)<1.05
            error('Nu must be greater than 1 when using GED errors');
        end
    elseif errorType==4 %GED
        if length(startingvals)~=(p+q+constant+1) || size(startingvals,2)~=1
            error('startingvals must be a column vector with p+q+CONSTANT+1 elements');
        end
        if startingvals(p+q+CONSTANT)<2.1
            error('Nu must be greater than 2.1 when using Skew T errors');
        end
        if startingvals(p+q+CONSTANT+1)<-.9 || startingvals(p+q+CONSTANT+1)>.9
            error('Lambda must be between -.9 and .9 when using Skew T errors');
        end
    end
    
    if any(startingvals(1:p+constant+q-1)<=0)
        error('All startingvals for omega, alpha and beta must be strictly greater than zero');
    end
    alphas = startingvals(constant+1:p+constant);
    betas = startingvals(constant+p+1:q+p+constant-1);
    if isempty(betas)
        betas = 0;
    end
    if (sum(alphas)+sum(betas))>=1
        error('The sum of the arch and garch coefficients must be less than 1');
    end
else
    startingvals=[];
end

%%%%%%%%%%%%%%%
% options
%%%%%%%%%%%%%%%
if nargin>7 && ~isempty(options)
    try
        optimset(options);
    catch OE
        error('options is not a valid minimization option structure');
    end
else
    %Setup the options in case of none provided
    options  =  optimset('fminunc');
    options  =  optimset(options , 'TolFun'      , 1e-005);
    options  =  optimset(options , 'TolX'        , 1e-005);
    options  =  optimset(options , 'Display'     , 'iter');
    options  =  optimset(options , 'Diagnostics' , 'on');
    options  =  optimset(options , 'LargeScale'  , 'off');
    options  =  optimset(options , 'MaxFunEvals' , 200*(2+p+q));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
