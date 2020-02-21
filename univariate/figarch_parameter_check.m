function [p,q,errorType,truncLag,startingvals,options]=figarch_parameter_check(epsilon,p,q,errorType,truncLag,startingvals,options)
% FIGARCH(Q,D,P) input validation.  Ensures that input parameters are conformable to what is
% expected.  
%
% USAGE:
%   [P,Q,ERRORTYPE,TRUNCLAG,STARTINGVALS,OPTIONS] = ...
%        figarch_parameter_check(EPSILON,P,Q,ERRORTYPE,TRUNCLAG,STARTINGVALS,OPTIONS);
%
% INPUTS:
%   See FIGARCH.
%
% OUTPUTS:
%   See FIGARCH.
%
% COMMENTS:
%
%  See also FIGARCH, FIGARCH_LIKELIHOOD, FIGARCH_STARTING_VALUES,
%  FIGARCH_TRANSFORM, FIGARCH_ITRANSFORM, FIGARCH_WEIGHTS 

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
% p
%%%%%%%%%%%%%%%
if isempty(p)
    p = 0;
end
if (length(p) > 1) || ~ismember(p,[0 1])
    error('P must be either 0 or 1.')
end

%%%%%%%%%%%%%%%
% q
%%%%%%%%%%%%%%%
if isempty(q)
    q = 0;
end
if (length(q) > 1) || ~ismember(q,[0 1])
    error('Q must be either 0 or 1.')
end

%%%%%%%%%%%%%%%
% errorType
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
% truncLag
%%%%%%%%%%%%%%%
if nargin<5 || isempty(truncLag)
    truncLag = 1000;
end
if truncLag<10 || floor(truncLag)~=truncLag || length(truncLag)>1
    error('TRUNCLAG must be a positive integer larger than 10.');
end



%%%%%%%%%%%%%%%%%
% starting values
%%%%%%%%%%%%%%%%%
if nargin>5 && ~isempty(startingvals)
    if size(startingvals,2)>size(startingvals,1)
        startingvals = startingvals';
    end
    %Validate starting vals, different if normal than if T or GED
    d = startingvals(2+p);
    if p
        phi = startingvals(2);
    end
    if q
        beta = startingvals(3+p);
    end
    omega = startingvals(1);
    if omega<=0
        error('Omega must be a positive scalar in STARTINGVALS.')
    end
    if d<=0 || d>=1
        error('d must be strictly between 0 and 1 in STARTINGVALS.')
    end
    if p && (phi>=((1-d)/2) || phi<=0)
        error('phi must satisfy 0<phi<(1-d)/2 with strict inequalities in STARTINGVALS.')
    end
    if q
        if p && (phi + d - beta)<=0 || beta<=0
            error('beta must satisfy 0<beta<phi + d with strict inequalities in STARTINGVALS.')
        elseif (d - beta)<=0 || beta<=0
            error('beta must satisfy 0<beta< d with strict inequalities in STARTINGVALS.')
        end
    end
    
    switch errorType
        case 1
            if length(startingvals)~=(p+q+2) || size(startingvals,2)~=1
                error('startingvals must be a column vector with 2+p+q elements');
            end
        case 2 %T
            if length(startingvals)~=(p+q+3) || size(startingvals,2)~=1
                error('startingvals must be a column vector with 2+p+q+1 elements');
            end
            if startingvals(p+q+3)<2.1
                error('Nu must be greater than 2.1 when using Students-T errors');
            end
        case 3 %GED
            if length(startingvals)~=(p+q+3) || size(startingvals,2)~=1
                error('startingvals must be a column vector with p+q+3 elements');
            end
            if startingvals(p+q+3)<1.05
                error('Nu must be greater than 1 when using GED errors');
            end
        case 4 %GED
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
        error('OPTIONS is not a valid minimization option structure');
    end
else
    %Setup the options in case of none provided
    options  =  optimset('fminunc');
    options  =  optimset(options , 'TolFun'      , 1e-005);
    options  =  optimset(options , 'TolX'        , 1e-005);
    options  =  optimset(options , 'Display'     , 'iter');
    options  =  optimset(options , 'Diagnostics' , 'on');
    options  =  optimset(options , 'LargeScale'  , 'off');
    options  =  optimset(options , 'MaxFunEvals' , 400*(2+p+q));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%