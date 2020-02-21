function options = realized_options(realizedFunction)
% Returns an option structure for use with realized_kernel and realized_variance_optimal_sampling containing the default values.
%
% USAGE:
%   [OPTIONS] = realized_options(REALIZEDFUNCTION)
%
% INPUTS:
%   REALIZEDFUNCTION - String containing the name of estimator being used.  Currently supported:
%                      'Kernel' - Realized kernel estimation (realized_kernel)
%                      'Multivariate Kernel' - Realized kernel estimation (realized_multivariate_kernel)
%                      'Optimal Sampling' - Realized variance using Bandi-Russell optimal sampling
%                        (realized_variance_optimal_sampling)
%                      'Twoscale' - Two-scale Realized variance 
%                      'Multiscale' - Multisclae Realized variance
%                      'QMLE' - QMLE estimation of QV
%                      'Preaveraging' - Preaveraged realized variance and related estimators
%                      
% OUTPUT:
%   OPTIONS - A structure that may include the fields (default values in parentheses)
%               kernel - One of the supported kernels ('nonflatparzen')
%                 Non-flat-top, weakly positive, n^(1/5) rate:
%                   'nonflatparzen' [RECOMMENDED] Parzen kernel applied to non-flat-top case
%                   'qs' Quadratic Spectral
%                   'fejer' Fejer kernel
%                   'thinf' Tukey-Hanning kernel with infinite lag order
%                   'bnhls' Kernel proposed by Barndorf-Neilsen, Hansen, Lunde and Shephard
%                 Flat-top, n^(1/4) rate:
%                   'parzen' Parzen's kernel
%                   'th1' Tukey-Hanning kernel with power 1
%                   'th2' Tukey-Hanning kernel with power 2
%                   'th5' Tukey-Hanning kernel with power 5
%                   'th16' Tukey-Hanning kernel with power 16
%                   'cubic','multiscale' 3rd order (cubic) kernel, which is asymptotically
%                     equivalent to multiscale realized variance
%                   '5thorder' 5th order Kernel
%                   '6thorder' 6th order Kernel
%                   '7thorder' 7th order Kernel
%                   '8thorder' 8th order Kernel
%                 Flat-top, n^(1/6) rate:
%                   'bartlett','twoscale' Bartlett kernel, which is asymptotically equivalent to the
%                     two-scale estimator
%                   '2ndorder' 2nd order (quadratic) kernel
%                   'epanechnikov' Epanechnikov kernel
%               endTreatment -  String value containing 'Jitter' or 'Stagger' ('Jitter')
%                 'Jitter'ing the end points means that the first q and the last q observations are
%                 pre-averaged before computing the kernel.  This is needed (theoretically) to
%                 remove end effects.  Once the end points are jittered the kernel is computed using
%                 all data to compute every realized autocovariance.  NOTE: 'Jitter' must be used if
%                 a non-flat top kernel is used.
%                 'Stagger' does not pre-average but uses a kernel computed from realized
%                 autocovariances computed from non-overlapping returns.
%               jitterLags - Number of lags to use if 'Jitter' is selected in endTreatment (2).  If
%               empty ([]) the optimal number of end points to jitter is selected automatically.
%               maxBandwidthPerc - Non-negative scalar between 0 and 1 maximum
%                 bandwidth to use as a percentage of the number of observations (.25)
%               maxBandwidth - Non-negative scalar maximum bandwidth to use as an absolut number ([])
%               bandwidth - Non-negative scale bandwidth to use when computing the realized kernel ([]).
%                 If empty the optimal bandwidth is estimated using a plug-in estimator.
%
%             The remaining fields are only used if bandwidth or jitterLags is estimated.
%             (i.e. OPTIONS.bandwidth = [] or OPTIONS.jitterLags = [])
%               useDebiasedNoise - Logical value indicating whether to use debiased noise or not
%                 when estimating the optimal bandwidth (false)
%               useAdjustedNoiseCount - Logical value indicating whether the count of the number of
%                 observations used in the noise estimate should include only non-zero returns (true)
%               medFrequencySamplingType - Sampling type to use when computing the realized
%                 variance and kernel used in the debiased noise estimate ('BusinessUniform').  See
%                 realized_price_filter for a description of the avilable sampling types
%               medFrequencySamplingInterval - Sampling interval to use when computing the realized
%                 variance and kernel used in the debiased noise estimate (390).  See
%                 realized_price_filter for the precise interpretation of sampling interval which
%                 depends on the selected sampling type.
%               medFrequencyKernel - Kernel to use when computing the debiased noise estimated ('parzen')
%               medFrequencyBandwidth - Bandwidth to use with medFrequencyKernel (5)
%               noiseVarianceSamplingType - Sampling type to use when computing the realized
%                 variance used in the noise estimator ('BusinessUniform').
%               noiseVarianceSamplingInterval - Sampling interval to use when computing the realized
%                 variance used in the noise estimator (120).
%               IQEstimationSamplingType -  Sampling type to use when computing the realized
%                 variance used in the IQ estimator ('BusinessUniform')
%               IQEstimationSamplingInterval - Sampling interval to use when computing the realized
%                 variance used in the IQ estimator (39)
%               theta - The scale to use in preaveraged estimators (1)
%
% COMMENTS:
%   These values have been calibrated to use with liquid NYSE TAQ data.  It may be necessary to
%   modify these values if computing realized kernels on illiquid data or on data which trades over
%   interval significantly different than 6.5 hours (such as 24-hour markets or where the market
%   is closed for part of the day)
%
%   Relevant options for each function:
%                                  |  Realized  |  Multivariate  |  Optimal   |  Twoscale,  |  QMLE  | Preaveraging
%                                  |  Kernel    |  Kernel        |  Sampling  |  Multiccale |        |      
%                                  |            |                |            |             |        |      
%   kernel                         |    x       |        x       |            |             |        |      
%   endTreatment                   |    x       |        x       |            |             |        |      
%   jitterLags                     |    x       |        x       |            |             |        |      
%   maxBandwidthPerc               |    x       |        x       |            |      x      |        |      
%   maxBandwidth                   |    x       |        x       |            |      x      |        |      
%   bandwidth                      |    x       |        x       |            |      x      |        |      
%   useDebiasedNoise               |    x       |        x       |     x      |      x      |        |      
%   useAdjustedNoiseCount          |    x       |        x       |     x      |      x      |        |      
%   medFrequencySamplingType       |    x       |        x       |     x      |      x      |    x   |     x
%   medFrequencySamplingInterval   |    x       |        x       |     x      |      x      |    x   |     x  
%   medFrequencyKernel             |    x       |        x       |     x      |      x      |        |     x      
%   medFrequencyBandwidth          |    x       |        x       |     x      |      x      |        |     x      
%   noiseVarianceSamplingType      |    x       |        x       |     x      |      x      |    x   |     x   
%   noiseVarianceSamplingInterval  |    x       |        x       |     x      |      x      |    x   |     x
%   IQEstimationSamplingType       |    x       |        x       |     x      |      x      |        |     x 
%   IQEstimationSamplingInterval   |    x       |        x       |     x      |      x      |        |     x 
%   theta                          |            |                |            |             |        |     x   
%
%  Note: Optimal Sampling uses default values of noiseVarianceSamplingType = 'BusinessTime' and
%        noiseVarianceSamplingInterval = 1
%
% EXAMPLE:
% % Generate the default options for estimating Realized Kernels
% options = realized_options('kernel')
%
%  See also REALIZED_KERNEL, REALIZED_VARIANCE_OPTIMAL_SAMPLING, REALIZED_TWOSCALE_VARIANCE,
%  REALIZED_MULTISCALE_VARIANCE, REALIZED_QMLE_VARIANCE


% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 5/1/2008

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>1
    error('0 or 1 input required.')
end
if nargin==0
    realizedFunction = 'kernel';
end
realizedFunction = lower(realizedFunction);
if ~ismember(realizedFunction,{'kernel','optimal sampling','twoscale','multiscale','multivariate kernel','qmle','preaveraging'})
    error('REALIZEDFUNCTION must be one of the listed types.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%


% Kernel
options.kernel = 'nonflatparzen';
% End Treatment
options.endTreatment = 'jitter';
% Max bandwidth percentage of observations
options.maxBandwidthPerc = .25;
% Max bandwidth (absolute)
options.maxBandwidth = [];
% Jitter Lags
options.jitterLags = 2;
% Bandwidth
options.bandwidth = [];
% These options only affect the estimation of bandwdith
% Use debiased noise
options.useDebiasedNoise = false;
% Use adjusted count
options.useAdjustedNoiseCount = true;
% Medium Frequency Sampling Type
options.medFrequencySamplingType = 'BusinessUniform';
% Medium Frequency Sampling Type
options.medFrequencySamplingInterval = 390;
% Medium Frequency Kernel
options.medFrequencyKernel = 'parzen';
% Medium Frequency Bandwidth
options.medFrequencyBandwidth = 5;
% Noise Variance Sampling Type
options.noiseVarianceSamplingType = 'BusinessUniform';
% Noise Variance Sampling Interval
options.noiseVarianceSamplingInterval = 120;
% IQ Estimation Sampling Type
options.IQEstimationSamplingType = 'BusinessUniform';
% IQ Estimation Sampling Interval
options.IQEstimationSamplingInterval = 39;
% Theta for preaveraging
options.theta = 1;


% Remove any unneeded fields
switch realizedFunction
    case {'kernel','multivariate kernel'}
        fieldList = {'kernel','endTreatment','jitterLags','maxBandwidthPerc','maxBandwidth',...
            'bandwidth','useDebiasedNoise','useAdjustedNoiseCount','medFrequencySamplingType',...
            'medFrequencySamplingInterval','medFrequencyKernel','medFrequencyBandwidth',...
            'noiseVarianceSamplingType','noiseVarianceSamplingInterval','IQEstimationSamplingType',...
            'IQEstimationSamplingInterval'};
    case {'optimal sampling','twoscale','multiscale'}
        fieldList = {'bandwidth','useDebiasedNoise','useAdjustedNoiseCount','medFrequencySamplingType',...
            'medFrequencySamplingInterval','medFrequencyKernel','medFrequencyBandwidth',...
            'noiseVarianceSamplingType','noiseVarianceSamplingInterval','IQEstimationSamplingType',...
            'IQEstimationSamplingInterval'};
        options.noiseVarianceSamplingType = 'BusinessTime';
        options.noiseVarianceSamplingInterval = 1;
        options.useAdjustedNoiseCount = false;
    case {'qmle'}
        fieldList = {'medFrequencySamplingType','medFrequencySamplingInterval',...
            'noiseVarianceSamplingType','noiseVarianceSamplingInterval'};
        options.noiseVarianceSamplingType = 'BusinessTime';
        options.noiseVarianceSamplingInterval = 1;
    case {'preaveraging'}
        fieldList = {'theta','medFrequencySamplingType','medFrequencySamplingInterval',...
            'medFrequencyKernel','medFrequencyBandwidth','noiseVarianceSamplingType',...
            'noiseVarianceSamplingInterval','IQEstimationSamplingType',...
            'IQEstimationSamplingInterval','useAdjustedNoiseCount'};
        options.noiseVarianceSamplingType = 'BusinessTime';
        options.noiseVarianceSamplingInterval = 1;
end


fields = fieldnames(options);
for i=1:length(fields)
    if ~ismember(fields(i),fieldList)
        options = rmfield(options,fields(i));
    end
end

