function [impulses,impulsesstd,hfig] = impulseresponse(y,constant,lags,leads,sqrttype,graph,het,uncorr)
% Computes impulse responses for a VAR(P) or irregular VAR(P) and standard
% errors under a variety of assumptions on the covariance of the errors:
%   * Conditionally Homoskedastic and Uncorrelated
%   * Conditionally Homoskedastic but Correlated
%   * Heteroskedastic but Conditionally Uncorrelated
%   * Heteroskedastic and Correlated
%
% USAGE:
%   [IMPULSES]=impulseresponse(Y,CONSTANT,LAGS,LEADS)
%   [IMPULSES,IMPULSESTD,HFIG]=impulseresponse(Y,CONSTANT,LAGS,LEADS,SQRTTYPE,GRAPH,HET,UNCORR)
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
%                        estimated error standard deviations.
%                     4 - Generalized impulse response of Peseran and Shin, which is equivalent to
%                        K reorderings of the covariance matrix where each variable is ordered
%                        first when computing the IR to the shock to that variable.
%                  If the input is a K by K positive definite matrix, it is used as the
%                  covariance square root for computing the impulse response function.
%   GRAPH       - [OPTIONAL] Logical variable (0 (no graph) or 1 (graph)) indicating whether the
%                   function should produce a bar plot of the sample autocorrelations and confidence
%                   intervals. Default is to produce a graphic (GRAPH=1).
%   HET         - [OPTIONAL] A scalar integer indicating the type of
%                 covariance estimator
%                    0 - Homoskedastic
%                    1 - Heteroskedastic [DEFAULT]
%   UNCORR      - [OPTIONAL] A scalar integer indicating the assumed structure of the error
%                   covariance matrix
%                    0 - Correlated errors  [DEFAULT]
%                    1 - Uncorrelated errors
%
% OUTPUTS:
%   IMPULSES    - K by K by LEADS matrix containing the impulse responses where IMPULSES(i,j,h)
%                   contains the impulse to Y(i) due to shock j at period h
%   IMPULSESSTD - K by K by LEADS matrix containing the impulse response where IMPULSES(i,j,h)
%                   contains the impulse to Y(i) due to shock j at period h
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
        het=[];
        uncorr=[];
    case 5
        graph=[];
        het=[];
        uncorr=[];
    case 6
        
        het=[];
        uncorr=[];
        
    case 7
        uncorr=[];
    case 8
        % Nothing
    otherwise
        error('4 to 8 inputs required')
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
    error('GRAPH must be a scalar, either 1 or 0.')
end
% Check het
if isempty(het)
    het=1;
end
if ~isscalar(het)
    error('HET must be either 0 or 1')
end
if ~ismember(het,[0 1])
    error('HET must be either 0 or 1')
end
% Check uncorr
if isempty(uncorr)
    uncorr=0;
end
if ~isscalar(uncorr)
    error('UNCORR must be either 0 or 1')
end
if ~ismember(uncorr,[0 1])
    error('UNCORR must be either 0 or 1')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the output of vectorar to do most of the work
[parameters,stderr,tstat,pval,const,conststd,r2,errors,s2,paramvec,vcv] = vectorar(y,constant,lags,het,uncorr);

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


% Finally to compute te G matrices for G(1), G(2), ... G(leads)
% To compute these I need impuse_response(0) impuse_response(-1)
% impuse_response(-2) .... impuse_response(-P)
maxlag=max(lags);
nP=K*(K*maxlag+1);


%Initialize the impulses and thier std dev holder
impulses=zeros(K,K,leads+1);
% Loop over the leads and the ei terms, this time scaling by sig12
for j=1:leads+1
    for i=1:K
        ei=zeros(K,1);
        ei(i)=1;
        impulses(:,i,j)=Xi(:,:,j)*sig12*ei;
    end
end


G=zeros(K*K,nP,leads);
for s=1:leads
    for i=1:s
        % The +1 is due to the "0" impulse being in position 1
        LHS = impulses(:,:,i-1+1);
        RHS=zeros(K,K*maxlag+1);
        for j=1:maxlag
            % Should append XiOrig(:,:,s-i+1) if s-i+1>=1
            if s-i+1>=1
                matrixToAdd = impulses(:,:,s-i+1)';
            else
                matrixToAdd = zeros(K);
            end
            RHS(:,(j-1)*K+2:j*K+1)=matrixToAdd;
        end
        G(:,:,s) = G(:,:,s) + kron(LHS,RHS);
    end
end

%Build a matrix of trues and falses that represent the actual parameters
%(true) and false parameters (false)
TFparameters=false(K,K*maxlag+constant);
if constant
    TFparameters(:,1)=true(K,1);
end
for i=1:maxlag
    if ismember(i,lags)
        TFparameters(:,(i-1)*K+constant+1:i*K+constant)=true(K,K);
    end
end
TFparameters=TFparameters';
TFparameters=TFparameters(:);
pl=find(TFparameters);

% Create the new vcv with 0's in all rows and cols that are not in the
% model
index=1;
vcv2=zeros(length(TFparameters));
for i=1:length(pl)
    temp=vcv(index,:);
    vcv2(pl(i),pl)=temp;
    index=index+1;
end


% Finally use VCV2 to compute the VCV of each Xi
vecXivcv=zeros(K*K,K*K,leads+1);
for i=1:size(G,3)
    vecXivcv(:,:,i+1)=G(:,:,i)*vcv2*G(:,:,i)';
end



% The final step is to compute the impulses, which are scaled by the
% appropriate sqrt of the covariance matrix,  The impulses will be
% Xi*sig12*ei for i=1,2,...,K
% I then need the covariance of this expression, which can be written as
% eye(K)*Xi*sig12*ei, which has the vec kron((sig12*ei)',eye(K))vec(Xi),
% and I have the covariance of vecXi above, so the variance of Xi*sig12*ei
% is kron((sig12*ei)',eye(K))*vecXivcv*kron((sig12*ei)',eye(K))' and the K
% standard deviations are on the diagonal

%Initialize the impulses and thier std dev holder
impulsesstd=zeros(K,K,leads+1);
% Loop over the leads and the ei terms, this time scaling by sig12
for j=1:leads+1
    for i=1:K
        % Remember Xi(:,:,maxlag) is really the first
        stdErr = sqrt(diag(vecXivcv(:,:,j)));
        stdErr = reshape(stdErr,K,K)';
        impulsesstd(:,:,j)= stdErr;
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
            h=plot((0:leads),squeeze(impulses(i,j,:)),(0:leads),squeeze(impulses(i,j,:))+1.96*squeeze(impulsesstd(i,j,:)),(0:leads),squeeze(impulses(i,j,:))-1.96*squeeze(impulsesstd(i,j,:)),[0 leads],[0 0]);
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
