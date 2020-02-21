function weights=realized_kernel_weights(options)
% Computes weights for Realized Kernels
%
% USAGE:
%   [WEIGHTS] = realized_kernel_weights(OPTIONS)
%
% INPUTS:
%   OPTIONS - A realized kernel options structure.  See help realized_options for details
%
% OUTPUTS:
%   WEIGHTS - H by 1 vector of kernel weights where H depends on the bandwidth and kernel
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL. See Barndorf-Nielsen,
%   Hansen, Lunde and Shephard (2008a, 2008b) for details about realized kernels and
%   their properties.
%
%  See also REALIZED_KERNEL, REALIZED_KERNEL_BANDWIDTH, 
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=1
    error('One input required.')
end
 
% List of flat top kernels
flatTopKernelList = {'bartlett','twoscale','2ndorder','epanechnikov',...
    'cubic','multiscale','5thorder','6thorder','7thorder','8thorder','parzen',...
    'th1','th2','th5','th16'};
 
% List of non flat top kernels
nonFlatTopKernelList = {'nonflatparzen','qs','fejer','thinf','bnhls'};
 
% Combined kernel list
kernelList = [flatTopKernelList nonFlatTopKernelList];
 
if ~isfield(options,'kernel') || ~ismember(options.kernel,kernelList)
    error('KERNEL must be a field of OPTIONS and one of the listed types.')
end
 
if ~isfield(options,'bandwidth') || options.bandwidth<0
    error('BANDWIDTH must be a field of OPTIONS and a non-negative scalar.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Take relevant fields
kernel = options.kernel;
bandwidth = options.bandwidth;
 
% A big lookup table
if bandwidth>0
    if ismember(options.kernel,flatTopKernelList)
        % If it is flat top then the first lag is actually 0 which produces a weights of 1
        H = round(bandwidth);
        x = (1:H)';
        x = (x-1)/H;
    else
        H = bandwidth;
    end
 
    switch lower(kernel)
        case {'bartlett','twoscale'}
            weights = 1-x;
        case {'2ndorder'}
            weights = 1-2*x+x.^2;
        case {'epanechnikov'}
            weights = 1-x.^2;
        case {'cubic','multiscale'}
            weights = 1-2*x.^2+2*x.^3;
        case '5thorder'
            weights = 1-10*x.^3+15*x.^4-6*x.^5;
        case '6thorder'
            weights = 1-15*x.^4+24*x.^5-10*x.^6;
        case '7thorder'
            weights = 1-21*x.^5+35*x.^6-15*x.^7;
        case '8thorder'
            weights = 1-28*x.^6+48*x.^7-21*x.^8;
        case 'parzen'
            weights = (1-6*x.^2+6*x.^3).*(x>=0 & x<=1/2) + 2*(1-x).^3.*(x>1/2 & x<1);
        case 'th1'
            weights = sin(pi/2*(1-x)).^2;
        case 'th2'
            weights = sin(pi/2*(1-x).^2).^2;
        case 'th5'
            weights = sin(pi/2*(1-x).^5).^2;
        case 'th16'
            weights = sin(pi/2*(1-x).^16).^2;
        case 'nonflatparzen'
            x = (1:H)';
            x = x./(H+1);
            weights = (1-6*x.^2+6*x.^3).*(x>=0 & x<=1/2) + 2*(1-x).^3.*(x>1/2 & x<1);
        case 'qs'
            % Truncate at 30 * H where H is BW
            x = unique(round(1:30*H))';
            x = x./(H+1);
            weights = 3./x.^2 .* (sin(x)./x - cos(x));
        case 'fejer'
            % Truncate at 30 * H where H is BW
            x = unique(round(1:30*H))';
            x = x./(H+1);
            weights = (sin(x)./x).^2;
        case 'thinf'
            % Truncate at 4 * H where H is BW
            x = unique(round(1:4*H))';
            x = x./(H+1);
            weights = sin((pi/2).*exp(-x)).^2;
        case 'bnhls'
            % Truncate at 10 * H where H is BW
            x = unique(round(1:10*H))';
            x = x./(H+1);
            weights = (1+x).*exp(-x);
        otherwise
            error('OPTIONS.KERNEL must be one of the listed types.')
    end
else
    weights=[];
end