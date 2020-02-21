function [p,o,q,errorType,userDelta,noUserDelta,startingvals,options]=aparch_parameter_check(data, p, o, q, errorType, userDelta, startingvals, options)
% APARCH(P,O,Q) input validation.  Ensures that the input parameters are
% conformable to what is expected.
%
% USAGE:
%   [P,O,Q,ERROR_TYPE,STARTINGVALS,OPTIONS] = ...
%        aparch_parameter_check(DATA,P,O,Q,ERROR_TYPE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   See APARCH.
%
% OUTPUTS:
%   See APARCH.
%
% COMMENTS:
%   See also APARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


%%%%%%%%%%%%%%%
% data
%%%%%%%%%%%%%%%
if size(data,2) > 1 || length(data)==1
    error('data series must be a column vector.')
elseif isempty(data)
    error('data is empty.')
end

%%%%%%%%%%%%%%%
% q
%%%%%%%%%%%%%%%
if (length(q) > 1) || any(q < 0) || isempty(q)
    error('q must be a non-negative scalar.')
end

%%%%%%%%%%%%%%%
% o
%%%%%%%%%%%%%%%
if (length(o) > 1) || any(o < 0) || isempty(o)
    error('o must be a non-negative scalar.')
end

%%%%%%%%%%%%%%%
% p
%%%%%%%%%%%%%%%
if (length(p) > 1) || any(p <  1) || isempty(p)
    error('p must be positive scalar.')
end

%%%%%%%%%%%%%%%%%%%%%%
% p and o restrictions
%%%%%%%%%%%%%%%%%%%%%%
if o>p
    error('O must be less than or equal to P')
end

%%%%%%%%%%%%%%%
% errorType
%%%%%%%%%%%%%%%
if nargin<5
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

%%%%%%%%%%%%%%
% userDelta
%%%%%%%%%%%%%%
if nargin>5 && ~isempty(userDelta)
    if ~isscalar(userDelta) || ~isreal(userDelta) || userDelta>4 || userDelta<.3
        error('USERDELTA must be a scalar between 0.3 and 4')
    end
    noUserDelta = false;
else
    userDelta=[];
    noUserDelta = true;
end

%%%%%%%%%%%%%%%%%
% starting values
%%%%%%%%%%%%%%%%%
if nargin>6
    if ~isempty(startingvals)
        %Validate starting vals, different if normal than if T or GED, also
        %epends on whether a userDelta was provided
        if errorType==1
            if length(startingvals)~=(p+o+q+1+noUserDelta) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+o+q+1 elements');
            end
        elseif errorType==2 %T
            if length(startingvals)~=(p+o+q+2+noUserDelta) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+o+q+2 elements');
            end
            if startingvals(p+o+q+2+noUserDelta)<2.1
                error('Nu must be greater than 2.1 when using Students-T errors');
            end
        elseif errorType==3 %GED
            if length(startingvals)~=(p+o+q+2+noUserDelta) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+o+q+2 elements');
            end
            if startingvals(p+o+q+2+noUserDelta)<1.05
                error('Nu must be greater than 1 when using GED errors');
            end
        elseif errorType==4 %SkewT
            if length(startingvals)~=(p+o+q+3+noUserDelta) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+o+q+3 elements');
            end
            if startingvals(p+o+q+2+noUserDelta)<2.1
                error('Nu must be greater than 2.1 when using Skew T errors');
            end
            if startingvals(p+o+q+3+noUserDelta)<-.9 || startingvals(p+o+q+3+noUserDelta)>.9
                error('Lambda must be between -.9 and .9 when using Skew T errors');
            end
        end
        if any(startingvals(1:p+1)<=0)
            error('All startingvals for omega and alpha must be strictly greater than zero');
        end
        if any(startingvals(p+2:p+o+1)>1) || any(startingvals(p+2:p+o+1)<-1)
            error('All starting values for gamma must be between -1 and 1')
        end
        if (sum(startingvals(2:p+1)) + sum(startingvals(p+o+2:p+o+q+1)))>=1
            error('The sum of the arch, garch and 0.5*tarch coefficients must be less than 1');
        end
        if noUserDelta
            if startingvals(p+o+q+2)>4 || startingvals(p+o+q+2)<=.3
                error('delta must be between .3 and 4')
            end
        end
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
    catch
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
if isempty(q)
    q=0;
end
if isempty(o)
    o=0;
end
