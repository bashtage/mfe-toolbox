function [text,AIC,BIC]=igarch_display(parameters,ll,vcv,data,p,q,errorType,igarchType,constant)
% Display parameters, tstats, pvals, log-likelihood and AIC/BIC
% from estimates of a IGARCH(P,O,Q) produced using igarch
%
% USAGE:
%   [TEXT] = igarch_display(PARAMETERS,LL,VCV,DATA,P,Q)
%   [TEXT,AIC,BIC] = igarch_display(PARAMETERS,LL,VCV,DATA,P,Q,ERRORTYPE,IGARCHTYE,CONSTANT)
%
% INPUTS:
%   PARAMETERS    - A CONSTANT+p+q-1 column vector of parameters with
%                   [omega alpha(1) ... alpha(p) beta(1) ... beta(q-1) [nu lambda]]'.
%   LL            - The log likelihood at the optimum
%   VCV           - Non-robust standard errors (inverse Hessian)
%   DATA          - A column of mean zero data
%   P             - Positive, scalar integer representing the number of
%                   symmetric innovations
%   Q             - Non-negative, scalar integer representing the number
%                   of lags of conditional variance (0 for ARCH)
%   ERRORTYPE     - [OPTIONAL] The error distribution used, valid types are:
%                     'NORMAL'    - Gaussian Innovations [DEFAULT]
%                     'STUDENTST' - T distributed errors
%                     'GED'       - Generalized Error Distribution
%                     'SKEWT'     - Skewed T distribution
%   IGARCHTYPE     - [OPTIONAL] The type of variance process, either
%                     1 - Model evolves in absolute values
%                     2 - Model evolves in squares [DEFAULT]
%   CONSTANT      - [OPTIONAL] Logical value indicating whether model
%                     should include a constant.  Default is true
%                     (include)
%
% OUTPUTS:
%   TEXT          - Character matrix with the formatted parameters of the model
%   AIC           - Aikake Information Criteria computed from the LL
%   BIC           - Schwartz/Bayesian Information Criteria computed from the LL
%
%  See also IGARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005



%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 6
        igarchType=2;
        errorType='NORMAL';
        constant = 1;
    case 7
        igarchType=2;
        constant = 1;
    case 8
        constant = 1;
    case 9
        % Nothing
    otherwise
        error('5 to 8 inputs required')
end
if isempty(igarchType)
    igarchType=2;
end
if isempty(errorType)
    errorType='NORMAL';
end
if isempty(constant)
    constant = 1;
end
% parameters, N by 1, real
if any(~isreal(parameters)) || size(parameters,2)~=1
    error('PARAMETERS must be a column vector.')
end
% LL
if ~isscalar(ll) || ~isreal(ll)
    error('LL must be a scalar.')
end
% VCV
if size(vcv,2)~=size(vcv,1) || any(min(eig(vcv))<=0) || size(vcv,1)~=length(parameters)
    error('VCV must be a square positive definite matrix compatible with PARAMETERS.')
end
% data
if any(~isreal(data)) || size(data,2)~=1
    error('DATA must be a T by 1 column vector.')
end
% p
if ~isscalar(p) || p<1 || floor(p)~=p
    error('P must be a postitive scalar integer')
end
% q
if ~isscalar(q) || q<0 || floor(q)~=q
    error('Q must be a positive scalar integer')
end
% errorType
if ~ischar(errorType)
    errorType=[];
end
switch errorType
    case 'NORMAL'
        errorType=1;
        extraP=0;
    case 'STUDENTST'
        errorType=2;
        extraP=1;
    case 'GED'
        errorType=3;
        extraP=1;
    case 'SKEWT'
        errorType=4;
        extraP=2;
    otherwise
        error('ERRORTYPE is not one of the supported distributions')
end
% igarchType
if ~ismember(igarchType,[1 2])
    error('IGARCHTYPE must be either 1 or 2')
end
if length(parameters)~=(constant+p+q-1+extraP)
    error('Size of PARAMETERS is not compatible with input P, O, Q')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%



if igarchType==1
    model_name = ['IAVARCH(' num2str(p) ',' num2str(q) ')'];
else
    model_name = ['IGARCH(' num2str(p) ',' num2str(q) ')'];
end

parameters_text=num2str(parameters,'%4.4f');
stderr=sqrt(diag(vcv));
stderr_text=num2str(stderr,'%4.4f');
tstats=parameters./stderr;
tstats_text=num2str(tstats,'%4.4f');
pvals=2-2*normcdf(abs(tstats));
pvals_text=num2str(pvals,'%4.4f');

% Add in the extra beta
finalBetaPos = constant+p+q;
archParameters = parameters(constant+1:constant+p+q-1);
finalBeta = 1 - sum(archParameters);

n = length(parameters);

parameters_text=[parameters_text(1:finalBetaPos-1,:);
    num2str(finalBeta,'%4.4f');
    parameters_text(finalBetaPos+1:n,:)];
stderr_text=[stderr_text(1:finalBetaPos-1,:);
    repmat('-',1,size(stderr_text,2));
    stderr_text(finalBetaPos+1:n,:)];
tstats_text=[tstats_text(1:finalBetaPos-1,:);
    repmat('-',1,size(tstats_text,2));
    tstats_text(finalBetaPos+1:n,:)];
pvals_text=[pvals_text(1:finalBetaPos-1,:);
    repmat('-',1,size(pvals_text,2));
    pvals_text(finalBetaPos+1:n,:)];





T=length(data);
AIC = -ll/T+2*length(parameters)/T;
BIC = -ll/T+log(T)*length(parameters)/T;




% Format the output
text=[];
text{1,1}= ' ';
text{2,1}= ' ';
text{3,1}=repmat('-',1,50);
text{4,1}=model_name;
text{5,1}=repmat('-',1,50);
text{6,1}=' ';
text{7,1}=['Loglikelihood: ' sprintf('%1.2f',ll)];
text{8,1}=['AIC: ' sprintf('%1.4f',AIC)];
text{9,1}=['BIC: ' sprintf('%1.4f',BIC)];
text{10,1}=' ';


for i=1:size(text,1);
    disp(text{i,1})
end


% Format the parameter table, need to right align everything
K=size(parameters_text,1);
for i=1:K
    N=size(parameters_text,2);
    if any(parameters_text(i,:)==' ')
        numSpace=sum(parameters_text(i,:)==' ');
        parameters_text(i,:)=[repmat(' ',1,numSpace) parameters_text(i,1:N-numSpace)];
    end
    N=size(stderr_text,2);
    if any(stderr_text(i,:)==' ')
        numSpace=sum(stderr_text(i,:)==' ');
        stderr_text(i,:)=[repmat(' ',1,numSpace) stderr_text(i,1:N-numSpace)];
    end
    N=size(tstats_text,2);
    if any(tstats_text(i,:)==' ')
        numSpace=sum(tstats_text(i,:)==' ');
        tstats_text(i,:)=[repmat(' ',1,numSpace) tstats_text(i,1:N-numSpace)];
    end
    N=size(pvals_text,2);
    if any(pvals_text(i,:)==' ')
        numSpace=sum(pvals_text(i,:)==' ');
        pvals_text(i,:)=[repmat(' ',1,numSpace) pvals_text(i,1:N-numSpace)];
    end
end
% Append column labels
labels={' Parameters','   Std. Err.','     T-stat','      P-val'};
cols = {parameters_text,stderr_text,tstats_text,pvals_text};
for i=1:length(labels);
    text1=labels{i};
    text2=cols{i};
    maxcols=max(size(text1,2),size(text2,2));
    if size(text1,2)<maxcols
        labels{i} = [repmat(' ',1,maxcols-size(text1,2)) text1 ];
    end
    if size(text2,2)<maxcols
        cols{i} = [repmat(' ',size(text2,1),maxcols-size(text2,2)) text2 ];
    end
end

% Finally variable names
variable_names=cell(length(parameters),1);
variable_names{1}='omega';
index=2;
for i=1:p
    variable_names{index}=['alpha(' num2str(i) ')'];
    index=index+1;
end
for i=1:q
    variable_names{index}=['beta(' num2str(i) ')'];
    index=index+1;
end
if errorType>1
    variable_names{index}='nu';
    index=index+1;
end
if errorType==4
    variable_names{index}='lambda';
end


maxVarNameLength=max(cellfun('length',variable_names));
for i=1:length(variable_names)
    variable_names{i} = [repmat(' ',1,maxVarNameLength-length(variable_names{i})) variable_names{i}];
end
variable_names=cell2mat(variable_names(:));

% Output the parameters
outmat=strvcat(' ',variable_names); %#ok<*VCAT>
for i=1:4
    outmat=[outmat repmat('  ',size(outmat,1),1) strvcat(labels{i},cols{i})]; %#ok<AGROW>
end


text2=[];
for i=1:length(text)
    text2=strvcat(text2,text{i});
end
text=strvcat(text2,outmat);

disp(outmat)
