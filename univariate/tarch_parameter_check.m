function [p,o,q,error_type,tarch_type,startingvals,options]=tarch_parameter_check(data, p, o, q, error_type, tarch_type, startingvals, options)
% TARCH(P,O,Q) input validation.  Ensures that the input parameters are
% conformable to what is expected.
%
% USAGE:
%   [P,O,Q,ERROR_TYPE,TARCH_TYPE,STARTINGVALS,OPTIONS] = ...
%        tarch_parameter_check(DATA,P,O,Q,ERROR_TYPE,TARCH_TYPE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   See TARCH.
%
% OUTPUTS:
%   See TARCH.
%
% COMMENTS:
%   See also TARCH

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

%%%%%%%%%%%%%%%
% error_type
%%%%%%%%%%%%%%%
if nargin<5
    error_type='NORMAL';
end

if isempty(error_type)
    error_type='NORMAL';
end

if strcmp(error_type,'NORMAL')
    error_type = 1;
elseif strcmp(error_type,'STUDENTST')
    error_type = 2;
elseif strcmp(error_type,'GED')
    error_type = 3;
elseif strcmp(error_type,'SKEWT')
    error_type = 4;
else
    error('error_type must be a string and one of: ''NORMAL'', ''STUDENTST'', ''GED'' or ''SKEWT''.');
end

%%%%%%%%%%%%%%%
% tarch_type
%%%%%%%%%%%%%%%
if nargin>5
    if isempty(tarch_type)
        tarch_type=2;
    end
    if length(tarch_type)>1
        error('tarch_type must be a scalar');
    end

    if ~(tarch_type==2 || tarch_type==1)
        error('tarch_type must be either 1 or 2')
    end
else
    tarch_type=2;
end

%%%%%%%%%%%%%%%%%
% starting values
%%%%%%%%%%%%%%%%%
if nargin>6
    if ~isempty(startingvals)
        %Validate starting vals, different if normal than if T or GED
        if error_type==1
            if length(startingvals)~=(p+o+q+1) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+o+q+1 elements');
            end
        elseif error_type==2 %T
            if length(startingvals)~=(p+o+q+2) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+o+q+2 elements');
            end
            if startingvals(p+o+q+2)<2.1
                error('Nu must be greater than 2.1 when using Students-T errors');
            end
        elseif error_type==3 %GED
            if length(startingvals)~=(p+o+q+2) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+o+q+2 elements');
            end
            if startingvals(p+o+q+2)<1.05
                error('Nu must be greater than 1 when using GED errors');
            end
        elseif error_type==4 %GED
            if length(startingvals)~=(p+o+q+3) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+o+q+3 elements');
            end
            if startingvals(p+o+q+2)<2.1
                error('Nu must be greater than 2.1 when using Skew T errors');
            end
            if startingvals(p+o+q+3)<-.9 || startingvals(p+o+q+3)>.9
                error('Lambda must be between -.9 and .9 when using Skew T errors');
            end
        end
        if any(startingvals(1:p+1)<=0)
            error('All startingvals for omega and alpha must be strictly greater than zero');
        end
        if (sum(startingvals(2:p+1))+ 0.5*sum(startingvals(p+2:p+o+1)) + sum(startingvals(p+o+2:p+o+q+1)))>=1
            error('The sum of the arch, garch and 0.5*tarch coefficients must be less than 1');
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
