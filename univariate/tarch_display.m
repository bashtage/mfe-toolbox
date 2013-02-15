function [text,AIC,BIC]=tarch_display(parameters,ll,vcv,data,p,o,q,errorType,tarchType)
% Display parameters, tstats, pvals, log-likelihood and AIC/BIC
% from estimates of a TARCH(P,O,Q) produced using tarch
%
% USAGE:
%   [TEXT] = tarch_display(PARAMETERS,LL,VCV,DATA,P,O,Q)
%   [TEXT,AIC,BIC] = tarch_display(PARAMETERS,LL,VCV,DATA,P,O,Q,ERRORTYPE,TARCHTYE)
%
% INPUTS:
%   PARAMETERS    - A 1+p+o+q column vector of parameters with
%                   [omega alpha(1) ... alpha(p) gamma(1) ... gamma(o) beta(1) ... beta(q) [nu lambda]]'.
%   LL            - The log likelihood at the optimum
%   VCV           - Non-robust standard errors (inverse Hessian)
%   DATA          - A column of mean zero data
%   P             - Positive, scalar integer representing the number of
%                   symmetric innovations
%   O             - Non-negative scalar integer representing the number
%                   of asymmetric innovations (0 for symmetric processes)
%   Q             - Non-negative, scalar integer representing the number
%                   of lags of conditional variance (0 for ARCH)
%   ERRORTYPE     - [OPTIONAL] The error distribution used, valid types are:
%                     'NORMAL'    - Gaussian Innovations [DEFAULT]
%                     'STUDENTST' - T distributed errors
%                     'GED'       - Generalized Error Distribution
%                     'SKEWT'     - Skewed T distribution
%   TARCHTYPE     - [OPTIONAL] The type of variance process, either
%                     1 - Model evolves in absolute values
%                     2 - Model evolves in squares [DEFAULT]
%
% OUTPUTS:
%   TEXT          - Character matrix with the formatted parameters of the model
%   AIC           - Aikake Information Criteria computed from the LL
%   BIC           - Schwartz/Bayesian Information Criteria computed from the LL
%
%  See also TARCH

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005



%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 5
        o=0;
        q=0;
        tarchType=2;
        errorType='NORMAL';
    case 6
        q=0;
        tarchType=2;
        errorType='NORMAL';
    case 7
        tarchType=2;
        errorType='NORMAL';
    case 8
        tarchType=2;
    case 9
    otherwise
        error('5 to 9 inputs required')
end
if nargin==7
elseif nargin==8
    errorType='NORMAL';
end
if isempty(tarchType)
    tarchType=2;
end
if isempty(errorType)
    errorType='NORMAL';
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
if size(vcv,2)~=size(vcv,1) || any(diag(vcv)<=0) || size(vcv,1)~=length(parameters)
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
% o
if ~isscalar(o) || o<0 || floor(o)~=o
    error('O must be a non-negative scalar integer')
end
% q
if ~isscalar(q) || q<0 || floor(q)~=q
    error('Q must be a non-negative scalar integer')
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
% tarchType
if ~ismember(tarchType,[1 2])
    error('TARCHTYPE must be either 1 or 2')
end
if length(parameters)~=(1+p+o+q+extraP)
    error('Size of PARAMETERS is not compatible with input P, O, Q')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%



if tarchType==1
    model_name = ['TARCH(' num2str(p) ',' num2str(o) ',' num2str(q) ')'];
else
    if o>0
        model_name = ['GJR-GARCH(' num2str(p) ',' num2str(o) ',' num2str(q) ')'];
    elseif q>0
        model_name = ['GARCH(' num2str(p)  ',' num2str(q) ')'];
    else
        model_name = ['ARCH(' num2str(p)  ')'];
    end
end

parameters_text=num2str(parameters,'%4.4f');
stderr=sqrt(diag(vcv));
stderr_text=num2str(stderr,'%4.4f');
tstats=parameters./stderr;
tstats_text=num2str(tstats,'%4.4f');
pvals=2-2*normcdf(abs(tstats));
pvals_text=num2str(pvals,'%4.4f');

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
for i=1:o
    variable_names{index}=['gamma(' num2str(i) ')'];
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
outmat=strvcat(' ',variable_names);
for i=1:4
    outmat=[outmat repmat('  ',size(outmat,1),1) strvcat(labels{i},cols{i})];
end


text2=[];
for i=1:length(text)
    text2=strvcat(text2,text{i});
end
text=strvcat(text2,outmat);

disp(outmat)
