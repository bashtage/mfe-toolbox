function  jitterLags = realized_kernel_jitter_lag_length(noiseEstimate,iqEstimate,kernel,N)
% Computes the optimal amount of end point jitter given a kernel, noise estimate, integrated
% quarticity estimate and number of observations.
%
% USAGE:
%   [JITTERLAGS] = realized_kernel_weights(NOISEESTIMATE,IQESTIMATE,KERNEL,N)
%
% INPUTS:
%   NOISEESTIMATE - Estimated variance of the noise present in a high freuqency asset price (see
%                     realized_kernel_select_lag_length)
%   IQESTIMATE    - Estimated integrated quarticity.  For the purposes of estimating the number of
%                     points to jitter it is often reasonable to use the square of a low frequence
%                     IV measure.
%   KERNEL        - String containing one of the supported kernel types:
%                     Non-flat-top, weakly positive, n^(1/5) rate:
%                     - 'nonflatparzen' [RECOMMENDED] Parzen kernel applied to non-flat-top case
%                     - 'qs' Quadratic Spectral
%                     - 'fejer' Fejer kernel
%                     - 'thinf' Tukey-Hanning kernel with infinite lag order
%                     - 'bnhls' Kernel proposed by Barndorf-Neilsen, Hansen, Lunde and Shephard
%                     Flat-top, n^(1/4) rate:
%                     - 'parzen' Parzen's kernel
%                     - 'th1' Tukey Hanning kernel with power 1
%                     - 'th2' Tukey Hanning kernel with power 2
%                     - 'th5' Tukey Hanning kernel with power 5
%                     - 'th16' Tukey Hanning kernel with power 16
%                     - 'cubic','multiscale' 3rd order (cubic) kernel, which
%                         is asymptotically equivalent to the multiscale estimator
%                     - '5thorder' 5th order Kernel
%                     - '6thorder' 6th order Kernel
%                     - '7thorder' 7th order Kernel
%                     - '8thorder' 8th order Kernel
%                     Flat-top, n^(1/6) rate:
%                     - 'bartlett','twoscale' Bartlett kernel, which is
%                         asymptotically equivalent to the two-scale estimator
%                     - '2ndorder' 2nd order (quadratic) kernel
%                     - 'epanechnikov' Epanechnikov kernel
%   N                 - Number of observation of the price
%
% OUTPUTS:
%   JITTERLAGS    - Estimated optimal number of lags to jitter the end points of a kernel.  If 1
%                     then no pre-averaging is needindicated
%
% COMMENTS:
%   This is a helper function for REALIZED_KERNEL. See Barndorf-Nielsen,
%   Hansen, Lunde and Shephard (2008a, 2008b) for details about the kernel and
%   their properties.
%
%  See also REALIZED_KERNEL_SELECT_LAG_LENGTH, REALIZED_VARIANCE, REALIZED_KERNEL, 
%  REALIZED_PRICE_FILTER, REALIZED_KERNEL_WEIGHTS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=4
    error('Four inputs required.')
end
if length(noiseEstimate)>1 || noiseEstimate<0
    error('NOISEESTIMATE must be a non-negative scalar.');
end
if length(iqEstimate)>1 || iqEstimate<0
    error('IQESTIMATE must be a non-negative scalar.');
end
if ~ismember(kernel,{'bartlett','twoscale','2ndorder','epanechnikov','cubic','multiscale','5thorder','6thorder','7thorder','8thorder','parzen','th1','th2','th5','th16','nonflatparzen','qs','fejer','thinf','bnhls'})
    error('KERNEL must be one of the listed types.')
end

if ~isscalar(N) || N<1 || floor(N)~=N
    error('N must be an integer greater than 1.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Look up table
switch lower(kernel)
    case {'bartlett','twoscale'}
        cStar = 2.28;
        k00 = 1/3;
        kernelType = 2;
    case '2ndorder'
        cStar = 3.42;
        k00 = 1/5;
        kernelType = 2;
    case 'epanechnikov'
        cStar = 2.46;
        k00 = 8/15;
        kernelType = 2;
    case {'cubic','multiscale'}
        cStar = 3.68;
        k00 = 0.371;
        k11 = 1.20;
        k22 = 12.0;
        kernelType = 1;
    case '5thorder'
        k00 = 0.391;
        k11 = 1.42;
        k22 = 17.1;
        cStar = 3.70;
        kernelType = 1;
    case '6thorder'
        k00 = 0.471;
        k11 = 1.55;
        k22 = 22.8;
        cStar = 3.97;
        kernelType = 1;
    case '7thorder'
        k00 = 0.533;
        k11 = 1.71;
        k22 = 31.8;
        cStar = 4.11;
        kernelType = 1;
    case '8thorder'
        k00 = 0.582;
        k11 = 1.87;
        k22 = 43.8;
        cStar = 4.31;
        kernelType = 1;
    case 'parzen'
        k00 = 0.269;
        k11 = 1.50;
        k22 = 24.0;
        cStar = 4.77;
        kernelType = 1;
    case 'th1'
        k00 = 0.375;
        k11 = 1.23;
        k22 = 12.1;
        cStar = 3.70;
        kernelType = 1;
    case 'th2'
        k00 = 0.219;
        k11 = 1.71;
        k22 = 41.7;
        cStar = 5.74;
        kernelType = 1;
    case 'th5'
        k00 = 0.097;
        k11 = 3.50;
        k22 = 489.0;
        cStar = 8.07;
        kernelType = 1;
    case 'th16'
        k00 = 0.032;
        k11 = 10.26;
        k22 = 14374.0;
        cStar = 39.16;
        kernelType = 1;
    case 'nonflatparzen'
        cStar = ((12)^2/0.269)^(1/5);
        k00 = .269;
        kernelType = 3;
    case 'qs'
        cStar = ((1/5)^2/(3*pi/5))^(1/5);
        k00 = 3*pi/5;
        kernelType = 3;
    case 'fejer'
        cStar = ((2/3)^2/(pi/3))^(1/5);
        k00 = pi/3;
        kernelType = 3;
    case 'thinf'
        cStar = ((pi^2/2)^2/(.52))^(1/5);
        k00 = 0.52;
        kernelType = 3;
    case 'bnhls'
        cStar = (1^2/(5/4))^(1/5);
        k00 = 5/4;
        kernelType = 3;
    otherwise
        error('KERNEL must be one of the listed types.')
end

        
% Compute the quantities necessary
if kernelType ==1
    % n^1/4 flat top
    % Approximate at the constant variance solution
    d=k00*k22/k11^2;
    fd = sqrt(1+sqrt(1+3*d));
    g = sqrt(k00*k11)*(1/fd + fd);
    avar = 16/3 * g * noiseEstimate * iqEstimate^(3/4);
    power = 1/4;
elseif kernelType == 2
    % n^1/6 flat top
    avar = 6 * cStar * k00 * noiseEstimate^(4/3) * iqEstimate^(2/3);
    power = 1/6;
elseif kernelType == 3
    % n^1/5 non flat top
    avar = 5 * cStar * k00 * noiseEstimate^(4/5) * iqEstimate^(4/5);
    power = 1/5;
end

% Problem is to min 8*noise^2*jitterLag^(-2) + avar * (N-jitterLag)^(-2*power)
% jittering over m data points, N observations in total
MSE = inf*ones(N-1,1);
for m=1:(N-1)
    MSE(m) = 8 * noiseEstimate^2 * m^(-2) + avar * (N-m)^(-2*power);
end

% Find the minimum, which may be 1
[temp,jitterLags] = min(MSE);