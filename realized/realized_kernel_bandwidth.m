function bandwidth = realized_kernel_bandwidth(noiseVariance, IQEstimate, options)
% Estimation of the optimal bandwidth for use in Realized Kernels
%
% USAGE:
%   [BANDWIDTH] = realized_kernel_select_bandwidth(NOISEVARIANCE,IQESTIMATE,OPTIONS)
%
% INPUTS:
%   NOISEVARIANCE   - Estimate of noise variance
%   IQESTIMATE      - Estimate of integrated quarticity
%   OPTIONS         - Realized kernel options structure.  See help realized_options
%
% OUTPUTS:
%   BANDWIDTH           - Optimal bandwidth computed using plug-in estimates of unknown quantities
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL.  See Barndorf-Nielsen, Hansen, Lunde and Shephard
%   (2008) for details about the optimal selection of bandwidth
%
%  See also REALIZED_KERNEL, REALIZED_NOISE_ESTIMATE
 
% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=3
    error('Three inputs required.')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Retrieve Kernel Specific information
switch lower(options.kernel)
    case {'bartlett','twoscale'}
        cStar = 2.28;
        kernelType = 2;
    case '2ndorder'
        cStar = 3.42;
        kernelType = 2;
    case 'epanechnikov'
        cStar = 2.46;
        kernelType = 2;
    case {'cubic','multiscale'}
        cStar = 3.68;
        kernelType = 1;
    case '5thorder'
        cStar = 3.70;
        kernelType = 1;
    case '6thorder'
        cStar = 3.97;
        kernelType = 1;
    case '7thorder'
        cStar = 4.11;
        kernelType = 1;
    case '8thorder'
        cStar = 4.31;
        kernelType = 1;
    case 'parzen'
        cStar = 4.77;
        kernelType = 1;
    case 'th1'
        cStar = 3.70;
        kernelType = 1;
    case 'th2'
        cStar = 5.74;
        kernelType = 1;
    case 'th5'
        cStar = 8.07;
        kernelType = 1;
    case 'th16'
        cStar = 39.16;
        kernelType = 1;
    case 'nonflatparzen'
        cStar = ((12)^2/0.269)^(1/5);
        kernelType = 3;
    case 'qs'
        cStar = ((1/5)^2/(3*pi/5))^(1/5);
        kernelType = 3;
    case 'fejer'
        cStar = ((2/3)^2/(pi/3))^(1/5);
        kernelType = 3;
    case 'thinf'
        cStar = ((pi^2/2)^2/(.52))^(1/5);
        kernelType = 3;
    case 'bnhls'
        cStar = (1^2/(5/4))^(1/5);
        kernelType = 3;
    otherwise
        error('KERNEL must be one of the listed types.')
end
 
 
% Compute the bandwidth
xiSquared = noiseVariance/sqrt(IQEstimate);
if xiSquared>1
    warning('oxfordRealized:excessiveLags','IQ estimaed o be close to 0.  Please check results');
    xiSquared = 1;
end
% Get the size of the filtered data 
nMax = options.filteredN;
 
if kernelType == 1
    bandwidth = cStar * (xiSquared)^(1/2) * nMax^(1/2);
elseif kernelType ==2
    bandwidth = cStar * (xiSquared)^(2/3) * nMax^(2/3);
elseif kernelType ==3
    bandwidth = cStar * (xiSquared)^(2/5) * nMax^(3/5);
end
