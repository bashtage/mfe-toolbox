function [parameters, LL, ht, VCVrobust, VCV, scores, diagnostics] = figarch(epsilon, p, q, errorType, truncLag, startingvals, options)
% FIGARCH(Q,D,P) parameter estimation for P={0,1} and Q={0,1} with different error distributions:
% Normal, Students-T, Generalized Error Distribution, Skewed T
%
% USAGE:
%   [PARAMETERS] = figarch(EPSILON,P,Q)
%   [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS] 
%                                       = figarch(EPSILON,P,Q,ERRORTYPE,TRUNCLAG,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   EPSILON       - T by 1 Column vector of mean zero residuals
%   P             - 0 or 1 indicating whether the autoregressive term is present in the model (phi)
%   Q             - 0 or 1 indicating whether the moving average term is present in the model (beta)
%   ERRORTYPE     - [OPTIONAL] The error distribution used, valid types are:
%                     'NORMAL'    - Gaussian Innovations [DEFAULT]
%                     'STUDENTST' - T distributed errors
%                     'GED'       - Generalized Error Distribution
%                     'SKEWT'     - Skewed T distribution
%   TRUNCLAG      - [OPTIONAL] Number of weights to compute in ARCH(oo) representation. 
%                      Default is 1000.
%   STARTINGVALS  - [OPTIONAL] A (2+p+q), plus 1 for STUDENTST OR GED (nu), plus 2 for SKEWT
%                      (nu,lambda), vector of starting values. [omega phi d beta [nu lambda]]'.  If
%                      not provided, FIGARCH_STARTING_VALUES attempts to find reasonable values.  
%   OPTIONS       - [OPTIONAL] A user provided options structure. Default options are below.
%
% OUTPUTS:
%   PARAMETERS    - A 2+p+q column vector of parameters with [omega phi d beta [nu lambda]]'.
%   LL            - The log likelihood at the optimum
%   HT            - The estimated conditional variances
%   VCVROBUST     - Robust parameter covariance matrix
%   VCV           - Non-robust standard errors (inverse Hessian)
%   SCORES        - Matrix of scores (# of params by t)
%   DIAGNOSTICS   - Structure of optimization output information.  Useful to check for convergence
%                     problems .
% COMMENTS:
%   The following (generally wrong) constraints are used:
%    (1) omega > 0
%    (2) 0<= d <= 1
%    (3) 0 <= phi <= (1-d)/2
%    (3) 0 <= beta <= d + phi
%    (5) nu>2 of Students T and nu>1 for GED
%    (6) -.99<lambda<.99 for Skewed T
%
%    The conditional variance, h(t), of a FIGARCH(1,d,1) process is modeled as follows:
%
%    h(t) = omega + [1-beta L - phi L  (1-L)^d] epsilon(t)^2 + beta * h(t-1)
%
%    where L is the lag operator which is estimated using an ARCH(oo) representation,
%
%    h(t) = (1-beta)^(-1) * omega + sum(lambda(i) * epsilon(t-i)^2)
%
%    where lambda(i) is a function of the fractional differencing parameter, phi and beta.
%
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
% See also TARCH, APARCH, EGARCH, AGARCH, FIGARCH_LIKELIHOOD, FIGARCH_PARAMETER_CHECK,
% FIGARCH_WEIGHTS FIGARCH_STARTING_VALUES, FIGARCH_TRANSFORM, FIGARCH_ITRANSFORM 

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/13/2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 3
        [p,q,errorType,truncLag,startingvals,options]=figarch_parameter_check(epsilon,p,q);
    case 4
        [p,q,errorType,truncLag,startingvals,options]=figarch_parameter_check(epsilon,p,q,errorType);
    case 5
        [p,q,errorType,truncLag,startingvals,options]=figarch_parameter_check(epsilon,p,q,errorType,truncLag);
    case 6
        [p,q,errorType,truncLag,startingvals,options]=figarch_parameter_check(epsilon,p,q,errorType,truncLag,startingvals);
    case 7
        [p,q,errorType,truncLag,startingvals,options]=figarch_parameter_check(epsilon,p,q,errorType,truncLag,startingvals,options);
    otherwise
        error('Number of inputs must be between 3 and 7');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Back Cast Value
backCastLength = max(floor(length(epsilon)^(1/2)),1);
backCastWeights = .05*(.9.^(0:backCastLength ));
backCastWeights = backCastWeights/sum(backCastWeights);
backCast = backCastWeights*((epsilon(1:backCastLength+1)).^2);
if backCast==0
    backCast=cov(epsilon);
end
% Back casting squared returns
epsilon2Augmented = [zeros(truncLag,1);epsilon.^2];
epsilon2Augmented(1:truncLag) = backCast;


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
[startingvals,nu,lambda,LLs,orderedParameters]=figarch_starting_values(startingvals,epsilon,epsilon2Augmented,p,q,errorType,truncLag);
%Finally, initialize the starting values
startingvals = [startingvals; nu; lambda];
%Transform the starting vals
[figarchParamsTransformed,nuTransformed,lambdaTransformed]=figarch_transform(startingvals,p,q,errorType);
%Re-append nu, lambda
startingvalsTransformed = [figarchParamsTransformed; nuTransformed; lambdaTransformed];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the parameters. Note the 1 in the last argument to indicate it
% is a constrained optimization

%LL0 is used to make sure the log likelihood improves
LL0=figarch_likelihood(startingvalsTransformed,p,q,epsilon,epsilon2Augmented,truncLag,errorType,true);

%Parameter estimation
[parameters,LL,exitflag,output]=fminunc('figarch_likelihood',startingvalsTransformed,options,p,q,epsilon,epsilon2Augmented,truncLag,errorType,true);
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
    [parameters,LL,exitflag,output]=fminunc('figarch_likelihood',parameters,options,p,q,epsilon,epsilon2Augmented,truncLag,errorType,true);
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
    while exitflag<=0 && index<=size(orderedParameters,1);
        %This condition checks that we haven't converged
        %OR that the best objective is worse than best grid search
        %Sort the original grid search log likelihoods and parameters
        
        startingvals=[orderedParameters(index,:)' ; nu; lambda];
        %Transform the starting vals
        [figarchParamsTransformed,nuTransformed,lambdaTransformed]=figarch_transform(startingvals,p,q,errorType);
        %Reappend nu
        startingvalsTransformed = [figarchParamsTransformed; nuTransformed; lambdaTransformed];
        
        LL0=figarch_likelihood(startingvalsTransformed,p,q,epsilon,epsilon2Augmented,truncLag,errorType,true);
        options.HessUpdate='bfgs';
        %Try the second set of starting values
        [parameters,LL,exitflag,output]=fminunc('figarch_likelihood',startingvalsTransformed,options,p,q,epsilon,epsilon2Augmented,truncLag,errorType,true);
        if  exitflag<=0 && LL<LL0
            %Again, if the LL improved, try more iterations
            %Increase the max iterations and max fun evals
            %Also switch to steepest descent
            options.MaxIter=2*options.MaxIter;
            options.MaxFunEvals=2*options.MaxFunEvals;
            options.HessUpdate='steepdesc';
            % Estimate the parameters.
            [parameters,LL,exitflag,output]=fminunc('figarch_likelihood',parameters,options,p,q,epsilon,epsilon2Augmented,truncLag,errorType,true);
        end
        %Save the parameter estimates
        robustParameters(index,:)=parameters';
        robustLL(index)=LL;
        %Increment the index
        index=index+1;
    end
    if index == size(orderedParameters,1);
        %save the best LL and parameters and break
        warning('MFEToolbox:Convergence','Convergence not achieved.  Use results with caution');
        [LL,index]=min(robustLL);
        parameters=robustParameters(index,:)';
    end
end

%Transform the parameters from the real line to the restricted space
[parameters,nu,lambda] = figarch_itransform(parameters,p,q,errorType);
parameters=[parameters;nu;lambda];
%Compute the log likelihood if needed
if nargout>1
    [LL, likelihoods, ht]=figarch_likelihood(parameters,p,q,epsilon,epsilon2Augmented,truncLag,errorType);
    LL=-LL;
end

%Compute standard errors using RobustVCV if needed.
if nargout>3
    nw=0; %No newey west on scores
    [VCVrobust,A,B,scores,hess]=robustvcv('figarch_likelihood',parameters,nw,p,q,epsilon,epsilon2Augmented,truncLag,errorType);
    T = size(epsilon,1);
    VCV=hess^(-1)/T;
end

%Report diagnostics in case requested
diagnostics.EXITFLAG=exitflag;
diagnostics.ITERATIONS=output.iterations;
diagnostics.FUNCCOUNT=output.funcCount;
diagnostics.MESSAGE=output.message;
