function [parameters, errors, SEregression, diagnostics, VCVrobust, VCV]=heterogeneousar(y,constant,p,nw,spec)
% Heterogeneous Autoregression parameter estimation
%
% USAGE:
%  [PARAMETERS] = heterogeneousar(Y,CONSTANT,P)
%  [PARAMETERS, ERRORS, SEREGRESSION, DIAGNOSTICS, VCVROBUST, VCV] 
%                                        = heterogeneousar(Y,CONSTANT,P,NW,SPEC)
%
% INPUTS:
%   Y            - A column of data
%   CONSTANT     - Scalar variable: 1 to include a constant, 0 to exclude
%   P            - A column vector or a matrix.
%                    If a vector, should include the indices to use for the lag length, such as in
%                    the usual case for monthly volatility data P=[1; 5; 22]. This indicates that
%                    the 1st lag, average of the first 5 lags, and the average of the first 22 lags
%                    should be used in estimation.  NOTE: When using the vector format, P MUST BE A
%                    COLUMN VECTOR to avoid ambiguity with the matrix format.  If P is a matrix, the
%                    values indicate the start and end points of the averages.  The above vector can
%                    be equivalently expressed as P=[1 1;1 5;1 22].  The matrix notation allows for
%                    the possibility of skipping lags, for example P=[1 1; 5 5; 1 22]; would have
%                    the 1st lag, the 5th lag and the average of lags 1 to 22.  NOTE: When using the
%                    matrix format, P MUST be # Entries by 2.  
%   NW           - [OPTIONAL] Number of lags to use when computing the long-run variance of the
%                    scores in VCVROBUST.  Default is 0. 
%   SPEC         - [OPTIONAL] String value indicating which representation to use in parameter
%                    estimation.  May be:
%                     'STANDARD' - Usual representation with overlapping lags
%                     'MODIFIED' - Modified representation with non-overlapping lags
%
% OUTPUTS:
%   PARAMETERS   - A 1+length(p) column vector of parameters with
%                  [constant har(1) ... har(P)]'
%   ERRORS       - A T by 1 length vector of errors from the regression with 0s in first max(max(P))
%                    places
%   SEREGRESSION - The standard error of the regressions
%   DIAGNOSTICS  - A structure of diagnostic information containing:
%                     P          - List of HAR lags used in estimation
%                     C          - Indicator if constant was included
%                     AIC        - Akaike Information Criteria for the estimated model
%                     SBIC       - Bayesian (Schwartz) Information Criteria for the
%                                   estimated model
%                     T          - Number of observations
%                     ADJT       - Length of sample used for estimation 
%                     ARROOTS    - The characteristic roots of the ARMA
%                                   process evaluated at the estimated parameters
%                     ABSARROOTS - The absolute value (complex modulus if
%                                    complex) of the ARROOTS
%   VCVROBUST    - Robust parameter covariance matrix, White if NW = 0,
%                    Newey-West if NW>0
%   VCV          - Non-robust standard errors (inverse Hessian)
%
% EXAMPLES:
%   Simulate data from a HAR model
%       y = armaxfilter_simulate(1000,1,22,[.1 .3/4*ones(1,4) .55/17*ones(1,17)])
%   Standard HAR with 1, 5 and 22 day lags
%       parameters = heterogeneousar(Y,1,[1 5 22]')
%   Standard HAR with 1, 5 and 22 days lags using matrix notation
%       parameters = heterogeneousar(Y,1,[1 1;1 5;1 22])
%   Standard HAR with 1, 5 and 22 day lags using the non-overlapping reparameterization
%       parameters = heterogeneousar(Y,1,[1 5 22]',[],'MODIFIED')
%   Standard HAR with 1, 5 and 22 day lags with Newey-West standard errors
%       [parameters, errors, seregression, diagnostics, vcvrobust, vcv] = ...
%                    heterogeneousar(Y,1,[1 5 22]',ceil(length(Y)^(1/3)))
%   Nonstandard HAR with lags 1, 2  and 10-22 day lags 
%       parameters = heterogeneousar(Y,1,[1 1;2 2;10 22])
%
% See also ARMAXFILTER, TARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/13/2009


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 3
        nw = 0;
        spec = 'STANDARD';
    case 4
        spec = 'STANDARD';
    case 5
        % Nothing
    otherwise
        error('3 to 5 inputs required');
end
%%%%%%%%%%%%%%%
% y
%%%%%%%%%%%%%%%
if size(y,2) > 1 || length(y)==1
    error('y series must be a column vector.')
elseif isempty(y)
    error('y is empty.')
end

%%%%%%%%%%%%%%%
% P
%%%%%%%%%%%%%%%
% P must be either a vector or a 2 by numP matrix, or possibly empty
if ~isempty(p)
    if size(p,1)>=1 && size(p,2)==2 % Then this will be interpreted as MATRIX format
        % Nothing to do now
    elseif min(size(p))==1 % Then this will be interpreted as MATRIX format
        p=[ones(size(p)) p];
    else
        error('P must be either a column vector or a # entries by 2 matrix.')
    end
end
% Eliminate any dupes
p = unique(p,'rows');
if min(min(p))<1 || any(any(p~=floor(p)))
    error('P must contain only positive integers.')
end


%%%%%%%%%%%%%%%
% Constant
%%%%%%%%%%%%%%%
if ~constant && isempty(p)
    error('At least one of CONSTANT or P must be nonempty.')
end
if ~ismember(constant,[0 1])
    error('CONSTANT must be 0 or 1.')
end

%%%%%%%%%%%%%%%
% NW
%%%%%%%%%%%%%%%
if isempty(nw)
    nw = 0;
end
if length(nw)>1 || floor(nw)~=nw || nw<0
    error('NW must be a non-negative integer.');
end

%%%%%%%%%%%%%%%
% Spec
%%%%%%%%%%%%%%%
if isempty(spec)
    spec = 'STANDARD';
end
spec = upper(spec);
if ~ismember(spec,{'STANDARD','MODIFIED'})
    error('SPEC must be either ''STANDARD'' or ''MODIFIED''.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Spec transform if needed
numP = size(p,1);
maxP = max(max(p));
if strcmp(spec,'MODIFIED')
    ind = int32(zeros(numP,maxP));
    for i=1:numP
        ind(i,p(i,1):p(i,2)) = 1;
    end
    used = [];
    % Transform ind to row echelon form
    for i=1:numP % Loop across columns
        notUsed = setdiff(1:numP,used);
        rows = ind(notUsed,:);
        if size(rows,1)>1
            [entry,position] = max(max(abs(rows)));
        else
            [entry,position] = max(abs(rows));
        end
        if entry>0
            [temp,rowSelected] = max(abs(rows(:,position))); %#ok<ASGLU>
            originalPosition = notUsed(rowSelected);
            used = [used originalPosition]; %#ok<AGROW>
            row = rows(rowSelected,:);
            leadignElement = row(position);
            row = row/leadignElement;
            % Find leading term location and Normalize to be 1
            for j=[1:originalPosition-1 originalPosition+1:numP]
                ind(j,:) = ind(j,:) - ind(j,position)*row;
                % Subtract i from j if same position leading term
            end
        end
    end
    % Finally clean up ind
    for i=1:numP
        element = ind(i,find(ind(i,:)~=0, 1 ));
        if ~isempty(element)
            ind(i,:) = ind(i,:)/element;
        end
    end
    %  Is amenable to transformation if max(sum(abs(ind)))==1
    newp = zeros(size(p));
    if max(sum(abs(ind)))==1
        for i=1:numP
            newp(i,:)=[find(ind(i,:), 1 ) find(ind(i,:), 1, 'last' )];
        end
        p = newp;
    else
        warning('MFEToolbox:IncompatibleInput','Input P is not compatible with the MODIFIED parameterization.  Using STANDARD parameterization');
    end
    p = sortrows(p);
end



[Y,X] = newlagmatrix(y,maxP,0);
T = length(Y);
newX = zeros(T,numP);
for i=1:numP
    newX(:,i) = mean(X(:,p(i,1):p(i,2)),2);
end
if constant
    newX = [ones(size(newX,1),1) newX];
end
numX = size(newX,2);
parameters = newX\Y;
errors = Y - newX*parameters;
SEregression = sqrt(errors'*errors/(T-numP));
A = newX'*newX/T;
Ainv = A\eye(size(A));
s = newX .* repmat(errors,1,numX);
B = covnw(s,nw,0);
VCVrobust = Ainv*B*Ainv/T;
VCV = SEregression.^2 * Ainv/T;
% Pad errors with 0s for return
errors = [zeros(length(y)-T,1);errors];



if nargout>=4
    diagnostics.P    = p;
    diagnostics.C    = constant;
    diagnostics.T    = length(y);
    diagnostics.adjT = T;
    diagnostics.spec = spec;
    [aic, sbic ]     = aicsbic(errors,constant,numX-constant,0);
    diagnostics.AIC  = aic;
    diagnostics.SBIC = sbic;
    
    if constant
        noConstParam = parameters(2:numX);
    else
        noConstParam = parameters;
    end
    coefficients = zeros(1,maxP);
    for i=1:numP
        coefficients(p(i,1):p(i,2)) = coefficients(p(i,1):p(i,2)) + noConstParam(i)/(p(i,2)-p(i,1)+1);
    end
    diagnostics.ARparameterization = coefficients;
    [arroots, absarroots]  = armaroots(coefficients,0,1:maxP,0);
    diagnostics.arroots    = arroots;
    diagnostics.absarroots = absarroots;
end
