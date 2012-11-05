function [error, errortext, sizeOut, varargout]=iscompatible(narg,varargin)
% Check whether a set of input parameters is compatible with a requested output size
%
% USAGE:
%   [ERROR, ERRORTEXT, SIZEOUT] = iscompatible(NARG,P1,P2,...,PK,S1,S2,...,SN)
%   [ERROR, ERRORTEXT, SIZEOUT] = iscompatible(NARG,P1,P2,...,PK,[S1 S2 ... SN])
%   [ERROR, ERRORTEXT, SIZEOUT, P1, P2, ... , PK] = iscompatible(NARG,P1,P2,...,PK,S1,S2,...,SN)
%   [ERROR, ERRORTEXT, SIZEOUT, P1, P2, ... , PK] = iscompatible(NARG,P1,P2,...,PK,[S1 S2 ... SN])
%
% INPUTS:
%   NARG      - Number of parameters to be checked, K above
%   Px        - Parameters to be chekced.  Must be as many as narg.  May be
%               either scalars of arrays
%   Sx        - Requested output size.  Either individual scalars or a vector
%
% OUTPUTS:
%   ERROR     - 1 if there is a problem with the input parameters, 0 otherwise
%   ERRORTEXT - A message detailing the nature of the problem
%   SIZEOUT   - The matched output size
%   Px        - Parameters transformed to all have the same dimension so that size(Px)=SIZEOUT
%
% COMMENTS:
%   Helper function for various distribution, probability, inverse
%   distribution and random number generators

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

% General Strategy
% 1. Sort the parametes form the size arguements
% 2. Find max dim of parameters, make sure all are either 1 by 1 or this
%    size; be careful if ndim is 2, then you want the max length
% 3. Figure out the output size, either a vector or a buch of scalars
% 4. Finally compare the parameter size to the SIZEOUT
% 5. Repmat the parameters is necessary

% Default error text
errortext='';
varargout=cell(narg);

% Error if there are not enough inputs
if length(varargin)<narg  || any(cellfun('isempty',varargin))
    sizeOut=[];
    error=1;
    errortext='Too few parameters.  All inputs must be nonempty.';
    return
end


% Parse the parameters
params=varargin(1:narg);
% Get parameter lengths
param_len=cellfun('length',params);

if all(param_len==1)
    % If they are all scalars, nothing to do
    param_size=[1 1];

elseif sum(param_len>1)==1
    % If only one is non-scalar, it's size becomes the common size
    param_size=size(params{param_len>1});

else
    % If more than one are scalar, get all teh sizes and make sure they are
    % the same
    param_non_scalar=find(param_len>1);
    params_non_scalar=params(param_non_scalar);
    for i=1:size(params_non_scalar)-1
        % This is an error if any two non-scalar parameters have different
        % sizes
        if ~isequal(size(params_non_scalar{i}),size(params_non_scalar{i+1}))
            sizeOut=[];
            error=1;
            errortext='Parameter size mismatch.  Must be either scalar or of a common size.';
            return
        end
    end
    %If they all have the same size, then the parameter size is just hte
    %size of the first non-scalar parameter
    param_size=size(params_non_scalar{1});
end


% Parse the requested output size
if length(varargin)>narg
    sizeOut=varargin(narg+1:length(varargin));
    if length(sizeOut)==1
        sizeOut=sizeOut{1};
        if ndims(sizeOut)==2 && length(sizeOut)==numel(sizeOut)
        else
            sizeOut=[];
            error=1;
            errortext=['Requested output size cannot be parsed.  Should be a ' ...
                'vector or series of scalars.'];
            return
        end
    else
        if all(cellfun('length',sizeOut)==1)
            sizeOut=cell2mat(sizeOut);
        else
            sizeOut=[];
            error=1;
            errortext=['Requested output size cannot be parsed.  Should be a ' ...
                'vector or series of scalars.'];
            return
        end
    end
else
    sizeOut=[];
end



% Finally, compare the parameter size to the requested output size
if prod(param_size)~=1 && ~isempty(sizeOut)
    if isequal(param_size,sizeOut)
        % If both are provided, are they equal
        error=0;
    else
        % If not, an error
        error=1;
        sizeOut=[];
        errortext='Requested output size and parameters are not of compatible sizes.';
        return
    end
elseif isempty(sizeOut)
    % If no sizeOut is requested, then the parameter size is the sizeOut
    sizeOut=param_size;
    error=0;
else
    % If all parameters are scalars, then and a sizeOut is provided, use that
    error=0;
end

% If the user wants transformed parameters to the common size, provide them
if nargout>3
    for i=1:narg;
        if length(varargin{i})==1
            varargout{i}=repmat(varargin{i},sizeOut);
        else
            varargout{i}=varargin{i};
        end
    end
end