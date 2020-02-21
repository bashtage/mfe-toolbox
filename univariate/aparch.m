function [parameters, LL, ht, VCVrobust, VCV, scores, diagnostics] = aparch(data, p, o, q, errorType, userDelta, startingvals, options)
% APARCH(P,O,Q) parameter estimation with different error distributions:
% Normal, Students-T, Generalized Error Distribution, Skewed T
%
% USAGE:
%   [PARAMETERS] = aparch(DATA,P,O,Q)
%   [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS] 
%                    = aparch(DATA,P,O,Q,ERRORTYPE,USERDELTA,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA         - A column of mean zero data
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   O            - Non-negative scalar integer representing the number of asymmetric innovations (0
%                    for symmetric processes) 
%   Q            - Non-negative, scalar integer representing the number of lags of conditional
%                    variance (0 for ARCH) 
%   ERRORTYPE    - [OPTIONAL] The error distribution used, valid types are:
%                    'NORMAL'    - Gaussian Innovations [DEFAULT]
%                    'STUDENTST' - T distributed errors
%                    'GED'       - Generalized Error Distribution
%                    'SKEWT'     - Skewed T distribution
%   USERDELTA    - [OPTIONAL] A scalar value between 0.3 and 4 to use for delta in the estimation.
%                    When the user provides a fixed value for delta, the vector of PARAMETERS has
%                    one less element.  This is useful for testing an unrestricted APARCH against
%                    TARCH or GJR-GARCH alternatives   
%   STARTINGVALS - [OPTIONAL] A (1+p+o+q+1), plus 1 for STUDENTST OR GED (nu),  plus 2 for SKEWT
%                    (nu,lambda), vector of starting values. 
%                  [omega alpha(1)...alpha(p) gamma(1)...gamma(o) beta(1)...beta(q) delta [nu lambda]]'.
%   OPTIONS      - [OPTIONAL] A user provided options structure. Default options are below.
%
% OUTPUTS:
%   PARAMETERS   - A 1+p+o+q+1 (+1 or 2) column vector of parameters with
%                  [omega alpha(1)...alpha(p) gamma(1)...gamma(o) beta(1)...beta(q) delta [nu lambda]]'.
%   LL           - The log likelihood at the optimum
%   HT           - The estimated conditional variances
%   VCVROBUST    - Robust parameter covariance matrix
%   VCV          - Non-robust standard errors (inverse Hessian)
%   SCORES       - Matrix of scores (# of params by t)
%   DIAGNOSTICS  - Structure of optimization output information.  Useful to check for convergence
%                    problems 
%
% COMMENTS:
%   The following (generally wrong) constraints are used:
%    (1) omega > 0
%    (2) alpha(i) >= 0 for i = 1,2,...,p
%    (3) 1<gamma<1 for i=1,...,o
%    (3) beta(i)  >= 0 for i = 1,2,...,q
%    (4) delta>.3
%    (5) sum(alpha(i) + beta(k)) < 1 for i = 1,2,...p and k=1,2,...,q
%    (6) nu>2 of Students T and nu>1 for GED
%    (7) -.99<lambda<.99 for Skewed T
%
%    The conditional variance, h(t), of a APARCH(P,O,Q) process is modeled as follows:
%
%    h(t)^(delta/2) = omega
%             + alpha(1)*(abs(r(t-1))+gamma(1)*r(t-1))^delta + ...
%               alpha(p)*(abs(r(t-p))+gamma(p)*r(t-p))^delta +
%               beta(1)*h(t-1)^(delta/2) +...+ beta(q)*h(t-q)^(delta/2)
%
%   Default Options
%     options  =  optimset('fmincon');
%     options  =  optimset(options , 'TolFun'      , 1e-005);
%     options  =  optimset(options , 'TolX'        , 1e-005);
%     options  =  optimset(options , 'Display'     , 'iter');
%     options  =  optimset(options , 'Diagnostics' , 'on');
%     options  =  optimset(options , 'LargeScale'  , 'off');
%     options  =  optimset(options , 'MaxFunEvals' , '400*numberOfVariables');
%
%  See also APARCH_LIKELIHOOD, APARCH_CORE, APARCH_PARAMETER_CHECK, APARCH_STARTING_VALUES,
%  APARCH_TRANSFORM, APARCH_ITRANSFORM 
%
%  You should use the MEX files (or compile if not using Win64 Matlab) as they provide speed ups of
%  approx 100 times relative to the m file 

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 4
        [p,o,q,errorType,userDelta,deltaIsEstimated,startingvals,options]=aparch_parameter_check(data, p, o, q);
    case 5
        [p,o,q,errorType,userDelta,deltaIsEstimated,startingvals,options]=aparch_parameter_check(data, p, o, q, errorType);
    case 6
        [p,o,q,errorType,userDelta,deltaIsEstimated,startingvals,options]=aparch_parameter_check(data, p, o, q, errorType, userDelta);        
    case 7
        [p,o,q,errorType,userDelta,deltaIsEstimated,startingvals,options]=aparch_parameter_check(data, p, o, q, errorType, userDelta,  startingvals);
    case 8
        [p,o,q,errorType,userDelta,deltaIsEstimated,startingvals,options]=aparch_parameter_check(data, p, o, q, errorType, userDelta,  startingvals, options);
    otherwise
        error('Number of inputs must be between 4 and 8');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial setup
m  =  max([p o q]);


% Augment the data with back casts to avoid costly memory allocations
% Because it is APARCH the backcase must be computed in every iteration
data_aug=[zeros(m,1);data];
abs_data_aug=[mean(abs(data))*ones(m,1);abs(data)];
%Compute the length of the augmented data
T     = size(data_aug,1);



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
%TARCH for starting values.
[startingvals,nu,lambda,LLs,ordered_parameters] = aparch_starting_values(startingvals,data,p,o,q,errorType,deltaIsEstimated);
%Finally, initialize the starting values
startingvals = [startingvals; nu; lambda];
%Transform the starting vals
[aparch_params_transformed,nu_transformed,lambda_transformed] = aparch_transform(startingvals,p,o,q,errorType,deltaIsEstimated);
%Re-append nu, lambda
startingvals_transformed = [aparch_params_transformed; nu_transformed; lambda_transformed];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the parameters. Note the 1 in the last argument to indicate it
% is a constrained optimization
%LL0 is used to make sure the log likelihood improves

[LL0,temp,ht0]=aparch_likelihood(startingvals_transformed,data_aug,abs_data_aug,p,o,q,errorType,T,deltaIsEstimated,userDelta,1);
%Parameter estimation
[parameters,LL,exitflag,output]=fminunc('aparch_likelihood',startingvals_transformed,options,data_aug,abs_data_aug,p,o,q,errorType,T,deltaIsEstimated,userDelta,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Estimation Robustness
% %This portion of the code is to make sure that the optimization converged
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %This is the case where the optimization did not converge, but improved on
% %the initial log likelihood
% if  exitflag<=0 && LL<LL0
%     %Try more iterations, only do more iterations if the final likelihood is
%     %actually better than the initial
%     
%     %Increase the max iterations and max fun evals
%     %Also switch to steepest descent
%     if ischar(options.MaxFunEvals)
%         options.MaxIter=2*100*length(parameters);
%     else
%         options.MaxIter=2*options.MaxIter;
%     end
%     if ischar(options.MaxFunEvals)
%         options.MaxFunEvals=4*100*length(parameters);
%     else
%         options.MaxFunEvals=2*options.MaxFunEvals;
%     end
%     options.HessUpdate='steepdesc';
%     % Estimate the parameters.
%     [parameters,LL,exitflag,output]=fminunc('aparch_likelihood',parameters,options,data_aug,abs_data_aug,p,o,q,errorType,back_cast,T,1);
% end
% 
% 
% 
% %If the optimization still hasn't converged, try other starting values
% if startingflag==0 && exitflag<=0
%     %Keep track of the final estimates, if nothing converges, we will
%     %return the 
%     robust_parameters(1,:)=parameters';
%     %Also keep the LL
%     robust_LL(1)=LL;
%     %Keep track of the iteration
%     index=2;
%     while exitflag<=0
%         %This condition checks that we haven't converged
%         %OR that the best objective is worse than best grid search
%         %Sort the original grid search log likelihoods and parameters
%         startingvals=ordered_parameters(index,:)';
%         startingvals=[startingvals;nu;lambda];
%         %Transform the starting vals
%         [garch_params_transformed,nu_transformed,lambda_transformed]=aparch_transform(startingvals,p,o,q,errorType);
%         %Reappend nu
%         startingvals_transformed = [garch_params_transformed; nu_transformed; lambda_transformed];
% 
%         LL0=aparch_likelihood(startingvals_transformed,data_aug,abs_data_aug,p,o,q,errorType,back_cast,T,1);
%         options.HessUpdate='bfgs';
%         %Try the second set of starting values
%         [parameters,LL,exitflag,output]=fminunc('aparch_likelihood',startingvals_transformed,options,data_aug,abs_data_aug,p,o,q,errorType,back_cast,T,1);
%         if  exitflag<=0 && LL<LL0
%             %Again, if the LL improved, try more iterations
%             %Increase the max iterations and max fun evals
%             %Also switch to steepest descent
%             options.MaxIter=2*options.MaxIter;
%             options.MaxFunEvals=2*options.MaxFunEvals;
%             options.HessUpdate='steepdesc';
%             % Estimate the parameters.
%             [parameters,LL,exitflag,output]=fminunc('aparch_likelihood',parameters,options,data_aug,abs_data_aug,p,o,q,errorType,back_cast,T,1);
%         end
%         %Save the parameter estimates
%         robust_parameters(index,:)=parameters';
%         robust_LL(index)=LL;
%         %Increment the index
%         index=index+1;
% 
%         if index>size(ordered_parameters,1);
%             %save the best LL and parameters and break
%             warning('Convergance not achieved.  Use results with caution');
%             [LL,index]=min(robust_LL);
%             parameters=robust_parameters(index,:)';
%             break
%         end
%     end
% end

%Transform the parameters from the real line to the restricted space
[parameters,nu,lambda]=aparch_itransform(parameters,p,o,q,errorType,deltaIsEstimated);
parameters=[parameters;nu;lambda];
%Compute the log likelihood if needed

fminunc('aparch_likelihood',startingvals_transformed,options,data_aug,abs_data_aug,p,o,q,errorType,T,deltaIsEstimated,userDelta,1);


if nargout>1
    [LL, likelihoods, ht]=aparch_likelihood(parameters,data_aug,abs_data_aug,p,o,q,errorType,T,deltaIsEstimated,userDelta);
    LL=-LL;
end

%Compute standard errors using RobustVCV if needed.
if nargout>3
    nw=0; %No newey west on scores
    [VCVrobust,A,B,scores,hess]=robustvcv('aparch_likelihood',parameters,nw,data_aug,abs_data_aug,p,o,q,errorType,T,deltaIsEstimated,deltaIsEstimated,userDelta);
    VCV=hess^(-1)/(T-m);
end

%Report diagnostics in case requested
diagnostics.EXITFLAG=exitflag;
diagnostics.ITERATIONS=output.iterations;
diagnostics.FUNCCOUNT=output.funcCount;
diagnostics.MESSAGE=output.message;
