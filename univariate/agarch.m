function [parameters, LL, ht, VCVrobust, VCV, scores, diagnostics] = agarch(epsilon, p, q, model_type, error_type, startingvals, options)
% AGARCH(P,Q) and NAGARCH(P,Q) with different error distributions:
% Normal, Students-T, Generalized Error Distribution, Skewed T
%
% USAGE:
%   [PARAMETERS] = agarch(EPSILON,P,Q)
%   [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS] 
%                                = agarch(EPSILON,P,Q,MODEL_TYPE,ERROR_TYPE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   EPSILON      - A column of mean zero data
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   Q            - Non-negative, scalar integer representing the number of lags of conditional
%                    variance (0 for ARCH-type model) 
%   MODEL_TYPE   - [OPTIONAL] The type of variance process, either
%                    'AGARCH'  - Asymmetric GARCH, Engle (1990) [DEFAULT]
%                    'NAGARCH' - Nonlinear Asymmetric GARCH, Engle & Ng (1993)
%   ERROR_TYPE   - [OPTIONAL] The error distribution used, valid types are:
%                    'NORMAL'    - Gaussian Innovations [DEFAULT]
%                    'STUDENTST' - T distributed errors
%                    'GED'       - Generalized Error Distribution
%                    'SKEWT'     - Skewed T distribution
%   STARTINGVALS - [OPTIONAL] A (2+p+q), plus 1 for STUDENTST OR GED (nu),  plus 2 for SKEWT
%                    (nu,lambda), vector of starting values. 
%                  [omega alpha(1) ... alpha(p) gamma beta(1) ... beta(q) [nu lambda]]'.
%   OPTIONS      - [OPTIONAL] A user provided options structure. Default options are below.
%
% OUTPUTS:
%   PARAMETERS   - A 2+p+q column vector of parameters with
%                  [omega alpha(1) ... alpha(p) gamma beta(1) ... beta(q) [nu lambda]]'.
%   LL           - The log likelihood at the optimum
%   HT           - The estimated conditional variances
%   VCVROBUST    - Robust parameter covariance matrix
%   VCV          - Non-robust standard errors (inverse Hessian)
%   SCORES       - Matrix of scores (# of params by t)
%   DIAGNOSTICS  - Structure of optimization output information.  Useful to check for convergence
%                     problems 
% COMMENTS:
%   The following (generally wrong) constraints are used:
%    (1) omega > 0
%    (2) alpha(i) >= 0 for i = 1,2,...,p
%    (3) beta(i)  >= 0 for i = 1,2,...,q
%    (4) -q(.01,EPSILON)<gamma<q(.99,EPSILON) for AGARCH 
%    (5) sum(alpha(i) + beta(k)) < 1 for i = 1,2,...p and k=1,2,...,q for
%    AGARCH and sum(alpha(i)*(1+gamma^2) + beta(k)) < 1 for NAGARCH
%    (6) nu>2 of Students T and nu>1 for GED
%    (7) -.99<lambda<.99 for Skewed T
%
%    The conditional variance, h(t), of a AGARCH(P,Q) process is given by:
%
%     h(t)  = omega
%             + alpha(1)*(r_{t-1}-gamma)^2 + ... + alpha(p)*(r_{t-p}-gamma)^2
%             + beta(1)*h(t-1) +...+ beta(q)*h(t-q)
%
%    The conditional variance, h(t), of a NAGARCH(P,Q) process is given by:
%
%     h(t)  = omega
%             + alpha(1)*(r_{t-1}-gamma*sqrt(h(t-1)))^2 + ... + alpha(p)*(r_{t-p}-gamma*sqrt(h(t-p)))^2 
%             + beta(1)*h(t-1) +...+ beta(q)*h(t-q)
%
%   Default Options
%     options  =  optimset('fminunc');
%     options  =  optimset(options , 'TolFun'      , 1e-005);
%     options  =  optimset(options , 'TolX'        , 1e-005);
%     options  =  optimset(options , 'Display'     , 'iter');
%     options  =  optimset(options , 'Diagnostics' , 'on');
%     options  =  optimset(options , 'LargeScale'  , 'off');
%     options  =  optimset(options , 'MaxFunEvals' , '200*numberOfVariables');
%
%  See also AGARCH_LIKELIHOOD, AGARCH_CORE, AGARCH_PARAMETER_CHECK, AGARCH_TRANSFORM, AGARCH_ITRANSFORM
%
%  You should use the MEX files (or compile if not using Win64 Matlab) as they provide speed ups of
%  approx 10 times relative to the m file 

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 7/12/2009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 3
        [p,q,model_type,error_type,startingvals,options]=agarch_parameter_check(epsilon, p, q);
    case 4
        [p,q,model_type,error_type,startingvals,options]=agarch_parameter_check(epsilon, p, q, model_type);
    case 5
        [p,q,model_type,error_type,startingvals,options]=agarch_parameter_check(epsilon, p, q, model_type, error_type);
    case 6
        [p,q,model_type,error_type,startingvals,options]=agarch_parameter_check(epsilon, p, q, model_type, error_type, startingvals);
    case 7
        [p,q,model_type,error_type,startingvals,options]=agarch_parameter_check(epsilon, p, q, model_type, error_type, startingvals, options);
    otherwise
        error('Number of inputs must be between 3 and 7');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initial setup
m  =  max([p q]);


% Augment the epsilon with local back casts to avoid costly memory allocations
back_cast_length = max(floor(length(epsilon)^(1/2)),1);
back_cast_weights = .05*(.9.^(0:back_cast_length ));
back_cast_weights = back_cast_weights/sum(back_cast_weights);
back_cast = back_cast_weights*((epsilon(1:back_cast_length+1)).^2);
if back_cast==0
     back_cast=cov(epsilon);
end
epsilon_augmented=[sqrt(back_cast)*ones(m,1);epsilon];
%Compute the length of the augmented epsilon
T = size(epsilon_augmented,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This flag is for the robustness check below, if the user supplies starting
%values, there will be no robustness check
[startingvals,nu,lambda]=agarch_starting_values(startingvals,epsilon,p,q,model_type,error_type);    
%Finally, initialize the starting values
startingvals = [startingvals; nu; lambda];
% Compute transform bounds
transform_bounds = quantile(epsilon,[.01 .99]);
%Transform the starting vals
[garch_params_transformed,nu_transformed,lambda_transformed]=agarch_transform(startingvals,p,q,model_type,error_type,transform_bounds);
%Re-append nu, lambda
startingvals_transformed = [garch_params_transformed; nu_transformed; lambda_transformed];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the parameters. Note the 1 in the last argument to indicate it
% is a constrained optimization

%LL0 is used to make sure the log likelihood improves
LL0=agarch_likelihood(startingvals_transformed,epsilon_augmented,p,q,model_type,error_type,transform_bounds ,back_cast,T,1);
%Parameter estimation
[parameters,LL,exitflag,output]=fminunc('agarch_likelihood',startingvals_transformed,options,epsilon_augmented,p,q,model_type,error_type,transform_bounds,back_cast,T,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Robustness
% This portion of the code is to make sure that the optimization converged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the case where the optimization did not converge, but improved on
% the initial log likelihood
if  exitflag<=0 && LL<LL0
    % Try more iterations, only do more iterations if the final likelihood is
    % actually better than the initial
    
    % Increase the max iterations and max fun evals
    % Also switch to steepest descent
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
    [parameters,LL,exitflag,output]=fminunc('agarch_likelihood',parameters,options,epsilon_augmented,p,q,model_type,error_type,transform_bounds,back_cast,T,1);
end

% Transform the parameters from the real line to the restricted space
[parameters,nu,lambda]=agarch_itransform(parameters,p,q,model_type,error_type,transform_bounds);
parameters=[parameters;nu;lambda];
% Compute the log likelihood if needed
if nargout>1
    [LL, likelihoods, ht]=agarch_likelihood(parameters,epsilon_augmented,p,q,model_type,error_type,transform_bounds,back_cast,T);
    LL=-LL;
end

%Compute standard errors using RobustVCV if needed.
if nargout>3
    nw=0; %No newey west on scores
    [VCVrobust,A,B,scores,hess]=robustvcv('agarch_likelihood',parameters,nw,epsilon_augmented,p,q,model_type,error_type,transform_bounds,back_cast,T);
    VCV=hess^(-1)/(T-m);
end

%Report diagnostics in case requested
diagnostics.EXITFLAG=exitflag;
diagnostics.ITERATIONS=output.iterations;
diagnostics.FUNCCOUNT=output.funcCount;
diagnostics.MESSAGE=output.message;
