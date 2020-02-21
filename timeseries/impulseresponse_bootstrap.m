function [impulses,impulseLowerCI,impulseUpperCI,hfig] = impulseresponse_bootstrap(y,constant,lags,leads,sqrttype,graph,bootstrap,B,w)
% Computes impulse responses for a VAR(P) or irregular VAR(P) and standard
% errors under a variety of assumptions on the covariance of the errors:
%   * Conditionally Homoskedastic and Uncorrelated
%   * Conditionally Homoskedastic but Correlated
%   * Heteroskedastic but Conditionally Uncorrelated
%   * Heteroskedastic and Correlated
%
% USAGE:
%   [IMPULSES]=impulseresponse(Y,CONSTANT,LAGS,LEADS)
%   [IMPULSES,LOWERCI,UPPERCI,HFIG]=impulseresponse(Y,CONSTANT,LAGS,LEADS,SQRTTYPE,GRAPH,BOOTSTRAP,B,W)
%
% INPUTS:
%   Y           - A T by K matrix of data
%   CONSTANT    - Scalar variable: 1 to include a constant, 0 to exclude
%   LAGS        - Non-negative integer vector representing the VAR orders to include in the model.
%   LEADS       - Number of leads to compute the impulse response function
%   SQRTTYPE    - [OPTIONAL] Either a scalar or a K by K positive definite matrix. This input
%                   determines the type of covariance decomposition used.  If it is a scalar if must
%                   be one  of:
%                     0 - Unit (unscaled) shocks, covariance assumed to be an identity matrix
%                     1 - [DEFAULT] Scaled but uncorrelated shocks. Scale is based on estimated
%                       error standard deviations.
%                     2 - Scaled and correlated shocks, Choleski decomposition. Scale is based on
%                       estimated error standard deviations.
%                     3 - Scaled and correlated shocks, spectral decomposition. Scale is based on
%                       estimated error standard deviations.
%                     4 - Generalized impulse response of Peseran and Shin, which is equivalent to
%                        K reorderings of the covariance matrix where each variable is ordered
%                        first when computing the IR to the shock to that variable.
%                  If the input is a K by K positive definite matrix, it is used as the
%                  covariance square root for computing the impulse response function.
%   GRAPH       - [OPTIONAL] Logical variable (0 (no graph) or 1 (graph)) indicating whether the
%                   function should produce a bar plot of the sample autocorrelations and confidence
%                   intervals. Default is to produce a graphic (GRAPH=1).
%  BOOTSTRAP    - [OPTIONAL] String value containing the type of bootstrap to use:
%                   'IID' for IID bootstrap (Default)
%                   'STATIONARY' for stationary bootstrap
%                   'BLOCK' for block bootstrap
%  B            - [OPTIONAL] Number of bootstrap replications to perform.  Default is 1000.
%  W            - [OPTIONAL] Positive integer containing the window length to use in the stationary
%                   or block bootstraps.  Ignored if BOOTSTRAP is 'IID'.  Default is 1.
%
% OUTPUTS:
%   IMPULSES    - K by K by LEADS matrix containing the impulse responses where IMPULSES(i,j,h)
%                   contains the impulse to Y(i) due to shock j at period h
%   LOWERCI     - K by K by LEADS matrix containing the boostrap lower bound where LOWERCI(i,j,h)
%                   contains the 2.5% bound for the impulse response of Y(i) due to shock j at period h
%   UPPERCI     - K by K by LEADS matrix containing the boostrap lower bound where LOWERCI(i,j,h)
%                   contains the 97.5% bound for the impulse response of Y(i) due to shock j at period h
%   HFIG        - Figure handle to the plot of the impulse response function
%
% COMMENTS:
%   Estimates a VAR including any lags.
%   y(:,t)' = CONST + P(1) * y(:,y-1) + P(2)*y(:,y-2) + ... + P(1)*y(:,t-K)'
%
%   where P(j) are K by K parameter matrices and CONST is a K by 1 parameter matrix (if CONSTANT==1)
%
% EXAMPLE:
%   To produce the IR for 12 leads form a VAR(1) with a constant
%       impulses = impulserepsonse(y,1,1,12)
%   To produce the IR for 12 leads form a VAR(3) without a constant
%       impulses  = impulserepsonse(y,0,[1:3],12)
%   To produce the IR for 12 leads form an irregular VAR(3) with only lags 1 and 3 with a constant
%       impulses  = impulserepsonse(y,1,[1:3],12)
%
% See also VECTORAR VECTORARVCV GRANGERCAUSE

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3.0    Date: 1/1/2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 4
        sqrttype=[];
        graph=[];
        bootstrap = [];
        B = [];
        w = [];
    case 5
        
        graph=[];
        bootstrap = [];
        B = [];
        w = [];
    case 6
        bootstrap = [];
        B = [];
        w = [];
    case 7
        B = [];
        w = [];
    case 8
        w = [];
    case 9
        % nothing
    otherwise
        error('4 to 9 inputs required')
end
% Check Y
if ndims(y)~=2
    error('Y must be T by K')
end
% Size of the cross-section
K=size(y,2);
% Check constant
if ~isscalar(constant)
    error('CONSTANT must be either 0 or 1')
end
if ~ismember(constant,[0 1])
    error('CONSTANT must be either 0 or 1')
end
% Check lags
if ndims(lags)~=2
    error('LAGS must be a vector of positive integers containing lags to include')
end
if size(lags,1)>size(lags,2)
    lags = lags';
end
if ~all(lags>0)
    error('LAGS must be a vector of positive integers containing lags to include')
end
if ~all(floor(lags)==lags)
    error('LAGS must be a vector of positive integers containing lags to include')
end
if length(lags)~=length(unique(lags))
    error('LAGS must be a vector of unique elements')
end
lags=sort(lags);
% Check leads
if ~isscalar(leads) || leads<1 || floor(leads)~=leads
    error('LEADS must be a positive scalar.')
end
% Check sqrttype
if isempty(sqrttype)
    sqrttype=1;
end
if isscalar(sqrttype)
    if ~ismember(sqrttype,[0 1 2 3 4])
        error('SQRTTYPE must be either a scalar (0,1,2,3,4) or a positive definite K by K matrix')
    end
    userCovSqrtProvided = false;
else
    if ndims(sqrttype)~=2
        error('SQRTTYPE must be either a scalar (0,1,2,3) or a positive definite K by K matrix')
    end
    if size(sqrttype,1)~=K || size(sqrttype,1)~=size(sqrttype,2)
        error('SQRTTYPE must be either a scalar (0,1,2,3) or a positive definite K by K matrix')
    end
    if min(eig(sqrttype))<=0
        error('SQRTTYPE must be either a scalar (0,1,2,3) or a positive definite K by K matrix')
    end
    userCovSqrtProvided = true;
end
% Check graph
if isempty(graph)
    graph=true;
end
if ~isscalar(graph) || ~ismember(graph,[0 1])
    error('GRAPH must be a scalar, either 1 or 0 (true or false).')
end
% check bootstrap
if isempty(bootstrap)
    bootstrap = 'iid';
end
bootstrap = lower(bootstrap);
if ~ismember(bootstrap,{'stationary','block','iid'})
    error('BOOTSTRAP must be either ''STATIONARY'', ''BLOCK'' or ''IID''');
end
% check B
if isempty(B)
    B = 1000;
end
if B<10 || floor(B)~=B
    error('B must be an integer greater than or equal to 10.  The Recommended value is 1000, and B should not be smaller than 100.')
end
% check w
if isempty(w)
    w = 1;
end
if w<1 || floor(w)~=w
    error('W must be a positive integer.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the output of vectorar to do most of the work
[parameters,stderr,tstat,pval,const,conststd,r2,errors,s2,paramvec,vcv] = vectorar(y,constant,lags);

%Need to know the longest lag
maxlag=max(lags);
% Find out which were not included in an irregular VAR
nolags=setdiff(1:maxlag,lags);

% Fill in any empty parameter matrices with 0s
for i=1:length(nolags)
    parameters{nolags(i)}=zeros(K);
end

if ~userCovSqrtProvided
    switch sqrttype
        case 0
            sig12 = eye(K);
        case 1
            sig12 = diag(sqrt(diag(s2)));
        case 2
            sig12 = chol(s2)';
        case 3
            sig12 = s2^(0.5);
        otherwise
            sig12 = zeros(K);
            for i=1:K
                order = [i setdiff(1:K,i)];
                s2Temp = s2(order,order);
                sTemp = chol(s2Temp)';
                [nothing,reOrder]=sort(order);
                sTemp = sTemp(reOrder,reOrder);
                sig12(:,i) = sTemp(:,i);
            end
    end
    
else
    sig12 = sqrttype;
end

% Set up the Xi matrics to hold the VMA(oo) representations
Xi=zeros(K,K,leads+1);
% Compute the VMA(oo) rep using recursion
for i=1:K
    this_impulse=zeros(leads+1,K);
    ei=zeros(K,1);
    ei(i)=1;
    shock=ei;
    this_impulse(1,:)=shock';
    for t=1:leads;
        for j=1:min(t,maxlag)
            this_impulse(t+1,:)=this_impulse(t+1,:)+(parameters{j}*this_impulse(t+1-j,:)')';
        end
    end
    Xi(:,i,:)=reshape(this_impulse',K,1,leads+1);
end


% Bootstrap the residuals and construct artificial time series
T = size(y,1);
switch bootstrap
    case 'stationary'
        indices = stationary_bootstrap((1:T-maxlag)',B,w);
    case 'block'
        indices = block_bootstrap((1:T-maxlag)',B,w);
    case 'iid'
        indices = ceil(rand(T-maxlag,B)*(T-maxlag));
end
% Fix indices to have T + maxlag
indices = [indices(T-3*maxlag+1:T-maxlag,:);indices];

initialValues = ceil(rand(B,1)*(T-(maxlag-1)))+(maxlag-1);
finalBootstrapImpulses = zeros(K,K,leads+1,B);
for i=1:B
    % FIX ME  Index problem
    bootstrapErrors = errors(indices(:,i),:);
    tempY = zeros(T+maxlag,K);
    tempY(1:maxlag,:) = y(initialValues(B)+(-maxlag+1:0),:);
    for t=maxlag+1:T+maxlag
        if constant
            tempY(t,:) = const;
        end
        for j=1:maxlag
            tempY(t,:) = tempY(t,:) + (parameters{j}*tempY(t-j,:)')';
        end
        tempY(t,:) = tempY(t,:) + bootstrapErrors(t-maxlag,:);
    end
    tempY = tempY(maxlag+1:T+maxlag,:);
    bootstrapParameters = vectorar(tempY,constant,lags);
    % Fill in any empty parameter matrices with 0s
    for ii=1:length(nolags)
        bootstrapParameters{nolags(ii)}=zeros(K);
    end
    % Set up the bootstrap Xi matrics to hold the VMA(oo) representations
    bootstrapXi=zeros(K,K,leads+1);
    % Compute the VMA(oo) rep using recursion
    for ii=1:K
        this_impulse=zeros(leads+1,K);
        ei=zeros(K,1);
        ei(ii)=1;
        shock=ei;
        this_impulse(1,:)=shock';
        for t=1:leads;
            for j=1:min(t,maxlag)
                this_impulse(t+1,:)=this_impulse(t+1,:)+(bootstrapParameters{j}*this_impulse(t+1-j,:)')';
            end
        end
        bootstrapXi(:,ii,:)=reshape(this_impulse',K,1,leads+1);
    end
    
    %Initialize the impulses and thier std dev holder
    bootstrapImpulses=zeros(K,K,leads+1);
    % Loop over the leads and the ei terms, this time scaling by sig12
    for j=1:leads+1
        for ii=1:K
            ei=zeros(K,1);
            ei(ii)=1;
            bootstrapImpulses(:,ii,j)=bootstrapXi(:,:,j)*sig12*ei;
        end
    end
    finalBootstrapImpulses(:,:,:,i) = bootstrapImpulses;
end
impulseLowerCI = quantile(finalBootstrapImpulses,.025,4);
impulseUpperCI = quantile(finalBootstrapImpulses,.975,4);

% The final step is to compute the impulses, which are scaled by the
% appropriate sqrt of the covariance matrix,  The impulses will be
% Xi*sig12*ei for i=1,2,...,K
% I then need the covariance of this expression, which can be written as
% eye(K)*Xi*sig12*ei, which has the vec kron((sig12*ei)',eye(K))vec(Xi),
% and I have the covariance of vecXi above, so the variance of Xi*sig12*ei
% is kron((sig12*ei)',eye(K))*vecXivcv*kron((sig12*ei)',eye(K))' and the K
% standard deviations are on the diagonal

%Initialize the impulses and thier std dev holder
impulses=zeros(K,K,leads+1);
% Loop over the leads and the ei terms, this time scaling by sig12
for j=1:leads+1
    for i=1:K
        ei=zeros(K,1);
        ei(i)=1;
        % Remember Xi(:,:,maxlag) is really the first
        impulses(:,i,j)=Xi(:,:,j)*sig12*ei;
    end
end

% Finally produce the plot if requested
if graph
    LB = zeros(K);
    UB = zeros(K);
    hfig=figure;
    set(hfig,'Position',[100 100 800 600])
    clf;
    for i=1:K
        for j=1:K
            subplot(K,K,(i-1)*K+j);
            h=plot((0:leads),squeeze(impulses(i,j,:)),(0:leads),squeeze(impulseLowerCI(i,j,:)),(0:leads),squeeze(impulseUpperCI(i,j,:)),[0 leads],[0 0]);
            set(h(1),'LineWidth',2,'Color',[.5 .5 1])
            set(h(2),'LineWidth',2,'Color',[.25 .25 .5],'LineStyle',':')
            set(h(3),'LineWidth',2,'Color',[.25 .25 .5],'LineStyle',':')
            set(h(4),'LineWidth',1,'Color',[0 0 0],'LineStyle','-')
            if i==1
                title(['e_' num2str(j)])
            end
            if j==1
                ylabel(['y_' num2str(i)])
            end
            axis tight
            AX=axis;
            spread=AX(4)-AX(3);
            AX(4)=AX(4)+.05*(spread);
            AX(3)=AX(3)-.05*(spread);
            LB(i,j) = AX(3);
            UB(i,j) = AX(4);
        end
    end
    UB = max(UB,[],2);
    LB = min(LB,[],2);
    for i=1:K
        for j=1:K
            subplot(K,K,(i-1)*K+j);
            AX = axis;
            AX(3) = LB(i);
            AX(4) = UB(i);
            axis(AX);
        end
    end
else
    hfig=[];
end
