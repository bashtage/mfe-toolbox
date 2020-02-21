function [parameters, LL, ht, VCVrobust, VCV, scores, diagnostics] = igarch(epsilon, p, q, errorType, igarchType, constant, startingvals, options)
% IGARCH(P,Q) parameter estimation with different error distributions
% Normal, Students-T, Generalized Error Distribution, Skewed T
% Estimation of IGARCH models  if IGARCHTYPE=2
% Estimation of IAVARCH if IGARCHTYPE=1
%
% USAGE:
%   [PARAMETERS] = igarch(EPSILON,P,Q)
%   [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS] 
%                         = igarch(EPSILON,P,Q,ERRORTYPE,IGARCHTYPE,CONSTANT,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   EPSILON      - A column of mean zero data
%   P            - Positive, scalar integer representing the number of innovations
%   Q            - Positive, scalar integer representing the number of lags of conditional variance 
%   ERRORTYPE    - [OPTIONAL] The error distribution used, valid types are:
%                    'NORMAL'    - Gaussian Innovations [DEFAULT]
%                    'STUDENTST' - T distributed errors
%                    'GED'       - Generalized Error Distribution
%                    'SKEWT'     - Skewed T distribution
%   IGARCHTYPE    - [OPTIONAL] The type of variance process, either
%                    1 - Model evolves in absolute values
%                    2 - Model evolves in squares [DEFAULT]
%   CONSTANT     - [OPTIONAL] Logical value indicating whether model should include a constant.
%                    Default is true (include). 
%   STARTINGVALS - [OPTIONAL] A P+Q, plus 1 for STUDENTST OR GED (nu), plus 2 for SKEWT
%                    (nu,lambda), vector of starting values. 
%                    [omega alpha(1) ... alpha(p) beta(1) ... beta(q) [nu lambda]]'.
%   OPTIONS      - [OPTIONAL] A user provided options structure. Default options are below.
%
% OUTPUTS:
%   PARAMETERS   - A P+Q column vector of parameters with
%                    [omega alpha(1) ... alpha(p) beta(1) ... beta(q-1) [nu lambda]]'.
%                    Note that the final beta is redundant and so excluded
%   LL           - The log likelihood at the optimum
%   HT           - The estimated conditional variances
%   VCVROBUST    - Robust parameter covariance matrix
%   VCV          - Non-robust standard errors (inverse Hessian)
%   SCORES       - Matrix of scores (# of params by t)
%   DIAGNOSTICS  - Structure of optimization output information.  Useful to check for convergence problems
%
% COMMENTS:
%   The parameters returned will always have P+Q elements, where the first
%   element is omega, and the remaining elements are
%   [alpha(1) alpha(2) ... alpha(P) beta(1) ... beta(Q-1)]
%   the value of the missing beta, beta(Q) is 1-sum(alpha) - sum(beta(1:Q-1))
%
%   The following (generally wrong) constraints are used:
%    (1) omega > 0 if CONSTANT
%    (2) alpha(i) >= 0 for i = 1,2,...,p
%    (3) beta(i)  >= 0 for i = 1,2,...,q
%    (4) sum(alpha(i) + beta(j)) = 1 for i = 1,2,...p and j = 1,2,...q
%    (5) nu>2 of Students T and nu>1 for GED
%    (6) -.99<lambda<.99 for Skewed T
%
%    The conditional variance, h(t), of an IGARCH(P,Q) process is modeled
%    as follows:
%
%     g(h(t)) = omega
%             + alpha(1)*f(r_{t-1}) + ... + alpha(p)*f(r_{t-p})+...
%             beta(1)*g(h(t-1)) +...+ beta(q)*g(h(t-q))
%
%     where f(x) = abs(x)  if IGARCHTYPE=1
%           g(x) = sqrt(x) if IGARCHTYPE=1
%           f(x) = x^2     if IGARCHTYPE=2
%           g(x) = x       if IGARCHTYPE=2
%
%   Default Options
%     options  =  optimset('fminunc');
%     options  =  optimset(options , 'TolFun'      , 1e-005);
%     options  =  optimset(options , 'TolX'        , 1e-005);
%     options  =  optimset(options , 'Display'     , 'iter');
%     options  =  optimset(options , 'Diagnostics' , 'on');
%     options  =  optimset(options , 'LargeScale'  , 'off');
%     options  =  optimset(options , 'MaxFunEvals' , '400*numberOfVariables');
%
% See also IGARCH_LIKELIHOOD, IGARCH_CORE, IGARCH_PARAMETER_CHECK, IGARCH_STARTING_VALUES,
% IGARCH_TRANSFORM, IGARCH_ITRANSFORM 
%
% You should use the MEX file for igarch_core (or compile if not using Win64 Matlab)
% as they provide speed ups of approx 100 times relative to the m file

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 3
        [p,q,errorType,igarchType,constant,startingvals,options]=igarch_parameter_check(epsilon, p, q);
    case 4
        [p,q,errorType,igarchType,constant,startingvals,options]=igarch_parameter_check(epsilon, p, q, errorType);
    case 5
        [p,q,errorType,igarchType,constant,startingvals,options]=igarch_parameter_check(epsilon, p, q, errorType, igarchType);
    case 6
        [p,q,errorType,igarchType,constant,startingvals,options]=igarch_parameter_check(epsilon, p, q, errorType, igarchType, constant);
    case 7
        [p,q,errorType,igarchType,constant,startingvals,options]=igarch_parameter_check(epsilon, p, q, errorType, igarchType, constant, startingvals);
    case 8
        [p,q,errorType,igarchType,constant,startingvals,options]=igarch_parameter_check(epsilon, p, q, errorType, igarchType, constant, startingvals, options);
    otherwise
        error('Number of inputs must be between 3 and 8');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initial setup
m  =  max([p q]);


%Augment the data with back casts to avoid costly memory allocations
if igarchType==1
    %fepsilon is f(epsilon), as above
    fepsilon   =  [mean(abs(epsilon))*ones(m,1) ; abs(epsilon)];
    
    % Local back casting
    backCastLength = max(floor(length(epsilon)^(1/2)),1);
    backCastWeights = .05*(.9.^(0:backCastLength ));
    backCastWeights = backCastWeights/sum(backCastWeights);
    backCast = backCastWeights*(abs(epsilon(1:backCastLength+1)));
    if backCast==0
        backCast=mean(abs(epsilon));
    end
else
    %fepsilon is f(epsilon), as above
    fepsilon   =  [cov(epsilon)*ones(m,1) ; epsilon.^2];

    % Local back casting
    backCastLength = max(floor(length(epsilon)^(1/2)),1);
    backCastWeights = .05*(.9.^(0:backCastLength ));
    backCastWeights = backCastWeights/sum(backCastWeights);
    backCast = backCastWeights*((epsilon(1:backCastLength+1)).^2);
    if backCast==0
        backCast=cov(epsilon);
    end
end
%Compute the length of the augmented epsilon
T     = size(fepsilon,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This flag is for the robustness check below, if the user supplies starting
%values, there will be no robustness check
if isempty(startingvals)
    startingflag=0;
else
    startingflag=1;
end
%Grid search for starting values.
[startingvals,nu,lambda,LLs,orderedParameters]=igarch_starting_values(startingvals,epsilon,fepsilon,p,q,T,errorType,igarchType,constant);
%Finally, initialize the starting values
startingvals = [startingvals; nu; lambda];
%Transform the starting vals
[garchParamsTransformed,nuTransformed,lambdaTransformed]=igarch_transform(startingvals,p,q,errorType,constant);
%Re-append nu, lambda
startingvalsTransformed = [garchParamsTransformed; nuTransformed; lambdaTransformed];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the parameters. Note the 1 in the last argument to indicate it
% is a constrained optimization

%LL0 is used to make sure the log likelihood improves
LL0=igarch_likelihood(startingvalsTransformed, epsilon, fepsilon, p, q, errorType, igarchType, constant, backCast, T, true);
%Parameter estimation
[parameters,LL,exitflag,output]=fminunc('igarch_likelihood',startingvalsTransformed,options, epsilon, fepsilon, p, q, errorType, igarchType, constant, backCast, T, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation Robustness
%This portion of the code is to make sure that the optimization converged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the case where the optimization did not converge, but improved on
%the initial log likelihood
if  exitflag<=0 && LL<LL0
    %Try more iterations, only do more iterations if the final likelihood is
    %actually better than the initial
    
    %Increase the max iterations and max fun evals
    %Also switch to steepest descent
    if ischar(options.MaxFunEvals)
        options.MaxIter=2*100*length(parameters);
    else
        options.MaxIter=2*options.MaxIter;
    end
    if ischar(options.MaxFunEvals)
        options.MaxFunEvals=4*100*length(parameters);
    else
        options.MaxFunEvals=2*options.MaxFunEvals;
    end
    options.HessUpdate='steepdesc';
    % Estimate the parameters.
    [parameters,LL,exitflag,output]=fminunc('igarch_likelihood',parameters,options, epsilon, fepsilon, p, q, errorType, igarchType, constant, backCast, T, true);
end



%If the optimization still hasn't converged, try other starting values
if startingflag==0 && exitflag<=0
    %Keep track of the final estimates, if nothing converges, we will
    %return the
    robustParameters(1,:)=parameters';
    %Also keep the LL
    robustLL = zeros(2,1);
    robustLL(1) = LL;
    %Keep track of the iteration
    index=2;
    while exitflag<=0
        %This condition checks that we haven't converged
        %OR that the best objective is worse than best grid search
        %Sort the original grid search log likelihoods and parameters
        
        startingvals=[orderedParameters(index,:)' ; nu; lambda];
        %Transform the starting vals
        [garchParamsTransformed,nuTransformed,lambdaTransformed]=igarch_transform(startingvals,p,q,errorType,constant);
        %Reappend nu
        startingvalsTransformed = [garchParamsTransformed; nuTransformed; lambdaTransformed];
        
        LL0=igarch_likelihood(startingvalsTransformed, epsilon, fepsilon, p, q, errorType, igarchType, constant, backCast, T, true);
        options.HessUpdate='bfgs';
        %Try the second set of starting values
        [parameters,LL,exitflag,output]=fminunc('igarch_likelihood',startingvalsTransformed,options, epsilon, fepsilon, p, q, errorType, igarchType, constant, backCast, T, true);
        if  exitflag<=0 && LL<LL0
            %Again, if the LL improved, try more iterations
            %Increase the max iterations and max fun evals
            %Also switch to steepest descent
            options.MaxIter=2*options.MaxIter;
            options.MaxFunEvals=2*options.MaxFunEvals;
            options.HessUpdate='steepdesc';
            % Estimate the parameters.
            [parameters,LL,exitflag,output]=fminunc('igarch_likelihood',parameters,options, epsilon, fepsilon, p, q, errorType, igarchType, constant, backCast, T, true);
        end
        %Save the parameter estimates
        robustParameters(index,:)=parameters';
        robustLL(index)=LL;
        %Increment the index
        index=index+1;
        
        if index>size(orderedParameters,1);
            %save the best LL and parameters and break
            warning('MFEToolbox:Convergence','Convergence not achieved.  Use results with caution');
            [LL,index]=min(robustLL);
            parameters=robustParameters(index,:)';
            break
        end
    end
end

%Transform the parameters from the real line to the restricted space
[parameters,nu,lambda]=igarch_itransform(parameters,p,q,errorType,constant);
parameters=[parameters;nu;lambda];
%Compute the log likelihood if needed
if nargout>1
    [LL, likelihoods, ht]=igarch_likelihood(parameters, epsilon, fepsilon, p, q, errorType, igarchType, constant, backCast, T);
    LL=-LL;
end

%Compute standard errors using RobustVCV if needed.
if nargout>3
    nw=0; %No newey west on scores
    [VCVrobust,A,B,scores,hess]=robustvcv('igarch_likelihood',parameters, nw, epsilon, fepsilon, p, q, errorType, igarchType, constant, backCast, T);
    VCV=hess^(-1)/(T-m);
end

%Report diagnostics in case requested
diagnostics.EXITFLAG=exitflag;
diagnostics.ITERATIONS=output.iterations;
diagnostics.FUNCCOUNT=output.funcCount;
diagnostics.MESSAGE=output.message;
