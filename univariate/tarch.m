function [parameters, LL, ht, VCVrobust, VCV, scores, diagnostics] = tarch(epsilon, p, o, q, error_type, tarch_type, startingvals, options)
% TARCH(P,O,Q) parameter estimation with different error distributions:
% Normal, Students-T, Generalized Error Distribution, Skewed T
% Estimation of ARCH or GARCH models if o=0 and tarch_type=2
% Estimation of TARCH or GJR asymmetric models if o>0 and tarch_type=1 or 2
%
% USAGE:
%   [PARAMETERS] = tarch(EPSILON,P,O,Q)
%   [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS] = 
%                                   tarch(EPSILON,P,O,Q,ERROR_TYPE,TARCH_TYPE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   EPSILON      - A column of mean zero data
%   P            - Positive, scalar integer representing the number of symmetric innovations
%   O            - Non-negative scalar integer representing the number of asymmetric innovations (0
%                    for symmetric processes)    
%   Q            - Non-negative, scalar integer representing the number of lags of conditional
%                    variance (0 for ARCH) 
%   ERROR_TYPE   - [OPTIONAL] The error distribution used, valid types are:
%                    'NORMAL'    - Gaussian Innovations [DEFAULT]
%                    'STUDENTST' - T distributed errors
%                    'GED'       - Generalized Error Distribution
%                    'SKEWT'     - Skewed T distribution
%   TARCH_TYPE   - [OPTIONAL] The type of variance process, either
%                    1 - Model evolves in absolute values
%                    2 - Model evolves in squares [DEFAULT]
%   STARTINGVALS - [OPTIONAL] A (1+p+o+q), plus 1 for STUDENTST OR GED (nu), plus 2 for SKEWT
%                    (nu,lambda), vector of starting values. 
%                     [omega alpha(1) ... alpha(p) gamma(1) ... gamma(o) beta(1) ... beta(q) [nu lambda]]'.
%   OPTIONS      - [OPTIONAL] A user provided options structure. Default options are below.
%
% OUTPUTS:
%   PARAMETERS   - A 1+p+o+q column vector of parameters with
%                  [omega alpha(1) ... alpha(p) gamma(1) ... gamma(o) beta(1) ... beta(q) [nu lambda]]'.
%   LL           - The log likelihood at the optimum
%   HT           - The estimated conditional variances
%   VCVROBUST    - Robust parameter covariance matrix
%   VCV          - Non-robust standard errors (inverse Hessian)
%   SCORES       - Matrix of scores (# of params by t)
%   DIAGNOSTICS  - Structure of optimization outputs and other values useful for functions calling TARCH.
% 
% COMMENTS:
% The following (generally wrong) constraints are used:
%   (1) omega > 0
%   (2) alpha(i) >= 0 for i = 1,2,...,p
%   (3) gamma(i) + alpha(i) > 0 for i=1,...,o
%   (3) beta(i)  >= 0 for i = 1,2,...,q
%   (4) sum(alpha(i) + 0.5*gamma(j) + beta(k)) < 1 for i = 1,2,...p and
%   j = 1,2,...o, k=1,2,...,q
%   (5) nu>2 of Students T and nu>1 for GED
%   (6) -.99<lambda<.99 for Skewed T
%
%   The conditional variance, h(t), of a TARCH(P,O,Q) process is modeled as follows:
%
%    g(h(t)) = omega
%            + alpha(1)*f(r_{t-1}) + ... + alpha(p)*f(r_{t-p})+...
%            + gamma(1)*I(t-1)*f(r_{t-1}) +...+ gamma(o)*I(t-o)*f(r_{t-o})+...
%            beta(1)*g(h(t-1)) +...+ beta(q)*g(h(t-q))
%
%     where f(x) = abs(x)  if tarch_type=1
%          g(x) = sqrt(x) if tarch_type=1
%          f(x) = x^2     if tarch_type=2
%          g(x) = x       if tarch_type=2
%
%   Default Options
%    options  =  optimset('fminunc');
%    options  =  optimset(options , 'TolFun'      , 1e-005);
%    options  =  optimset(options , 'TolX'        , 1e-005);
%    options  =  optimset(options , 'Display'     , 'iter');
%    options  =  optimset(options , 'Diagnostics' , 'on');
%    options  =  optimset(options , 'LargeScale'  , 'off');
%    options  =  optimset(options , 'MaxFunEvals' , '400*numberOfVariables');
%
%  See also TARCH_LIKELIHOOD, TARCH_CORE, TARCH_PARAMETER_CHECK, TARCH_STARTING_VALUES,
%  TARCH_TRANSFORM, TARCH_ITRANSFORM 
%
%  You should use the MEX files (or compile if not using Win64 Matlab) as they provide speed ups of
%  approx 100 times relative to the m file.

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 4
        [p,o,q,error_type,tarch_type,startingvals,options]=tarch_parameter_check(epsilon, p, o, q);
    case 5
        [p,o,q,error_type,tarch_type,startingvals,options]=tarch_parameter_check(epsilon, p, o, q, error_type);
    case 6
        [p,o,q,error_type,tarch_type,startingvals,options]=tarch_parameter_check(epsilon, p, o, q, error_type, tarch_type);
    case 7
        [p,o,q,error_type,tarch_type,startingvals,options]=tarch_parameter_check(epsilon, p, o, q, error_type, tarch_type, startingvals);
    case 8
        [p,o,q,error_type,tarch_type,startingvals,options]=tarch_parameter_check(epsilon, p, o, q, error_type, tarch_type, startingvals, options);
    otherwise
        error('Number of inputs must be between 4 and 8');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Initial setup
m  =  max([p o q]);


%Augment the data with back casts to avoid costly memory allocations
if tarch_type==1
    %fepsilon is f(epsilon), as above
    fepsilon   =  [mean(abs(epsilon))*ones(m,1) ; abs(epsilon)];
    %fIepsilon is fepsilon*(epsilon<0)
    fIepsilon   =  [0.5*mean(abs(epsilon))*ones(m,1) ; abs(epsilon).*(epsilon<0)];
    
    % Local back casting
    back_cast_length = max(floor(length(epsilon)^(1/2)),1);
    back_cast_weights = .05*(.9.^(0:back_cast_length ));
    back_cast_weights = back_cast_weights/sum(back_cast_weights);
    back_cast = back_cast_weights*(abs(epsilon(1:back_cast_length+1)));
    if back_cast==0
        back_cast=mean(abs(epsilon));
    end
else
    %fepsilon is f(epsilon), as above
    fepsilon   =  [mean(epsilon.^2)*ones(m,1) ; epsilon.^2];
    %fIepsilon is fepsilon*(epsilon<0)
    fIepsilon   =  [0.5*mean(epsilon.^2)*ones(m,1) ; epsilon.^2.*(epsilon<0)];
    % Local back casting
    back_cast_length = max(floor(length(epsilon)^(1/2)),1);
    back_cast_weights = .05*(.9.^(0:back_cast_length ));
    back_cast_weights = back_cast_weights/sum(back_cast_weights);
    back_cast = back_cast_weights*((epsilon(1:back_cast_length+1)).^2);
    if back_cast==0
        back_cast=mean(epsilon.^2);
    end
end
epsilon_augmented=[zeros(m,1);epsilon];
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
[startingvals,nu,lambda,~,ordered_parameters]=tarch_starting_values(startingvals,epsilon_augmented,fepsilon,fIepsilon,p,o,q,T,error_type,tarch_type);
%Finally, initialize the starting values
startingvals = [startingvals; nu; lambda];
%Transform the starting vals
[garch_params_transformed,nu_transformed,lambda_transformed]=tarch_transform(startingvals,p,o,q,error_type);
%Re-append nu, lambda
startingvals_transformed = [garch_params_transformed; nu_transformed; lambda_transformed];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the parameters. Note the 1 in the last argument to indicate it
% is a constrained optimization

%LL0 is used to make sure the log likelihood improves
LL0=tarch_likelihood(startingvals_transformed,epsilon_augmented,fepsilon,fIepsilon,p,o,q,error_type,tarch_type,back_cast,T,1);
%Parameter estimation
[parameters,LL,exitflag,output]=fminunc('tarch_likelihood',startingvals_transformed,options,epsilon_augmented,fepsilon,fIepsilon,p,o,q,error_type,tarch_type,back_cast,T,1);

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
    [parameters,LL,exitflag,output]=fminunc('tarch_likelihood',parameters,options,epsilon_augmented,fepsilon,fIepsilon,p,o,q,error_type,tarch_type,back_cast,T,1);
end



%If the optimization still hasn't converged, try other starting values
if startingflag==0 && exitflag<=0
    %Keep track of the final estimates, if nothing converges, we will
    %return the
    robust_parameters(1,:)=parameters';
    %Also keep the LL
    robust_LL = zeros(2,1);
    robust_LL(1) = LL;
    %Keep track of the iteration
    index=2;
    while exitflag<=0
        %This condition checks that we haven't converged
        %OR that the best objective is worse than best grid search
        %Sort the original grid search log likelihoods and parameters
        
        startingvals=[ordered_parameters(index,:)' ; nu; lambda];
        %Transform the starting vals
        [garch_params_transformed,nu_transformed,lambda_transformed]=tarch_transform(startingvals,p,o,q,error_type);
        %Reappend nu
        startingvals_transformed = [garch_params_transformed; nu_transformed; lambda_transformed];
        
        LL0=tarch_likelihood(startingvals_transformed,epsilon_augmented,fepsilon,fIepsilon,p,o,q,error_type,tarch_type,back_cast,T,1);
        options.HessUpdate='bfgs';
        %Try the second set of starting values
        [parameters,LL,exitflag,output]=fminunc('tarch_likelihood',startingvals_transformed,options,epsilon_augmented,fepsilon,fIepsilon,p,o,q,error_type,tarch_type,back_cast,T,1);
        if  exitflag<=0 && LL<LL0
            %Again, if the LL improved, try more iterations
            %Increase the max iterations and max fun evals
            %Also switch to steepest descent
            options.MaxIter=2*options.MaxIter;
            options.MaxFunEvals=2*options.MaxFunEvals;
            options.HessUpdate='steepdesc';
            % Estimate the parameters.
            [parameters,LL,exitflag,output]=fminunc('tarch_likelihood',parameters,options,epsilon_augmented,fepsilon,fIepsilon,p,o,q,error_type,tarch_type,back_cast,T,1);
        end
        %Save the parameter estimates
        robust_parameters(index,:)=parameters';
        robust_LL(index)=LL;
        %Increment the index
        index=index+1;
        
        if index>size(ordered_parameters,1);
            %save the best LL and parameters and break
            warning('MFEToolbox:Convergence','Convergence not achieved.  Use results with caution');
            [LL,index]=min(robust_LL);
            parameters=robust_parameters(index,:)';
            break
        end
    end
end

%Transform the parameters from the real line to the restricted space
[parameters,nu,lambda]=tarch_itransform(parameters,p,o,q,error_type);
parameters=[parameters;nu;lambda];
%Compute the log likelihood if needed
if nargout>1
    [LL, ~, ht]=tarch_likelihood(parameters,epsilon_augmented,fepsilon,fIepsilon,p,o,q,error_type,tarch_type,back_cast,T);
    LL=-LL;
end

%Compute standard errors using RobustVCV if needed.
if nargout>3
    nw=0; %No newey west on scores
    [VCVrobust,A,~,scores,hess]=robustvcv('tarch_likelihood',parameters,nw,epsilon_augmented,fepsilon,fIepsilon,p,o,q,error_type,tarch_type,back_cast,T);
    VCV=hess^(-1)/(T-m);
    diagnostics.A = A;
end

%Report diagnostics in case requested
diagnostics.EXITFLAG=exitflag;
diagnostics.ITERATIONS=output.iterations;
diagnostics.FUNCCOUNT=output.funcCount;
diagnostics.MESSAGE=output.message;
diagnostics.m = m;
diagnostics.T = T;
diagnostics.fdata = fepsilon;
diagnostics.fIdata = fIepsilon;
diagnostics.back_cast = back_cast;