function [parameters, LL, ht, VCVrobust, VCV, scores, diagnostics] = egarch(data, p, o, q, error_type, startingvals, options)
% EGARCH(P,O,Q) parameter estimation with different error distributions:
% Normal, Students-T, Generalized Error Distribution, Skewed T
%
% USAGE:
%   [PARAMETERS] = egarch(DATA,P,O,Q)
%   [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS] 
%                                 = egarch(DATA,P,O,Q,ERROR_TYPE,STARTINGVALS,OPTIONS)
%
% INPUTS:
%   DATA          - A column of mean zero data
%   P             - Positive, scalar integer representing the number of symmetric innovations
%   O             - Non-negative scalar integer representing the number of asymmetric innovations (0
%                     for symmetric processes) 
%   Q             - Non-negative, scalar integer representing the number of lags of conditional
%                   variance (0 for ARCH) 
%   ERROR_TYPE    - [OPTIONAL] The error distribution used, valid types are:
%                     'NORMAL'    - Gaussian Innovations [DEFAULT]
%                     'STUDENTST' - T distributed errors
%                     'GED'       - Generalized Error Distribution
%                     'SKEWT'     - Skewed T distribution
%   STARTINGVALS  - [OPTIONAL] A (1+p+o+q), plus 1 for STUDENTST OR GED (nu), 
%                   plus 2 for SKEWT (nu,lambda), vector of starting values.
%                   [omega alpha(1)...alpha(p) gamma(1)...gamma(o) beta(1)...beta(q) [nu lambda]]'.
%   OPTIONS       - [OPTIONAL] A user provided options structure. Default options are below.
%
% OUTPUTS:
%   PARAMETERS    - A 1+p+o+q column vector of parameters with
%                   [omega alpha(1)...alpha(p) gamma(1)...gamma(o) beta(1)...beta(q) [nu lambda]]'.
%   LL            - The log likelihood at the optimum
%   HT            - The estimated conditional variances
%   VCVROBUST     - Robust parameter covariance matrix
%   VCV           - Non-robust standard errors (inverse Hessian)
%   SCORES        - Matrix of scores (# of params by t)
%   DIAGNOSTICS   - Structure of optimization output information.  Useful to check for convergence
%                     problems 
%
% COMMENTS:
%   (1) Roots of the characteristic polynomial of beta are restricted to be less than 1
%
%   The conditional variance, h(t), of an EGARCH(P,O,Q) process is modeled as follows:
%
%   ln(h(t)) = omega
%            + alpha(1)*(abs(e_{t-1})-C) + ... + alpha(p)*(abs(e_{t-p})-C)+...
%            + gamma(1)*e_{t-1} +...+ e_{t-o} +...
%              beta(1)*ln(h(t-1)) +...+ beta(q)*ln(h(t-q))
%
%       where: ln is natural log
%              e_t = r_t/sqrt(h_t)
%              C   = 1/sqrt(pi/2)
%
%   Default Options
%    options  =  optimset('fmincon');
%    options  =  optimset(options , 'TolFun'      , 1e-005);
%    options  =  optimset(options , 'TolX'        , 1e-005);
%    options  =  optimset(options , 'Display'     , 'iter');
%    options  =  optimset(options , 'LargeScale'  , 'off');
%    options  =  optimset(options , 'MaxFunEvals' , 200*(2+p+q));
%    options  =  optimset(options , 'MaxSQPIter'  , 500);
%    options  =  optimset(options , 'Algorithm'   ,'active-set');
%
%  See also EGARCH_LIKELIHOOD, EGARCH_CORE, EGARCH_PARAMETER_CHECK, EGARCH_STARTING_VALUES, 
%           EGARCH_TRANSFORM, EGARCH_ITRANSFORM, EGARCH_NLCOM
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
        [p,o,q,error_type,startingvals,options]=egarch_parameter_check(data, p, o, q);
    case 5
        [p,o,q,error_type,startingvals,options]=egarch_parameter_check(data, p, o, q, error_type);
    case 6
        [p,o,q,error_type,startingvals,options]=egarch_parameter_check(data, p, o, q, error_type, startingvals);
    case 7
        [p,o,q,error_type,startingvals,options]=egarch_parameter_check(data, p, o, q, error_type, startingvals, options);
    otherwise
        error('Number of inputs must be between 4 and 8');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial setup
m  =  max([p o q]);

%Augment the data with backcasts to avoid costly memory allocations
back_cast_length = max(floor(length(data)^(1/2)),1);
back_cast_weights = .05*(.9.^(0:back_cast_length ));
back_cast_weights = back_cast_weights/sum(back_cast_weights);
back_cast = back_cast_weights*(data(1:back_cast_length+1).^2);
if back_cast==0
    back_cast=log(cov(data));
else
    back_cast=log(back_cast);
end
% absdata   =  [mean(abs(data))*ones(m,1) ; abs(data)];
data_augmented=[zeros(m,1);data];

%Compute the length of the augmented data
T     = size(data_augmented,1);
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
%Search for starting values.
[startingvals,nu,lambda,LLs,ordered_parameters]=egarch_starting_values(startingvals,data_augmented,p,o,q,T,error_type);
%Finally, initialize the starting values
startingvals = [startingvals; nu; lambda];
%Transform the starting vals
startingvals=egarch_transform(startingvals,p,o,q,error_type);
%LL0 is used to make sure the log-likelihood improves
LL0=egarch_likelihood(startingvals,data_augmented,p,o,q,error_type,back_cast,T,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up some constraints to prevent the estimation from going awry.
LB=ones(size(startingvals))*-inf;
%Need to keep alpha and gamma from making the process crazy
LB(2:p+o+1)=-3;
UB=ones(size(startingvals))*inf;
UB(2:p+o+1)=3;


% Estimate the parameters. Note the 1 in the last argument to indicate it
% is a constrained optimization
[parameters,LL,exitflag,output]=fmincon('egarch_likelihood',startingvals,[],[],[],[],LB,UB,'egarch_nlcon',options,data_augmented,p,o,q,error_type,back_cast,T,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation Robustness
%This portion of the code is to make sure that the optimization converged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the case where the optimization did not converge, but improved on
%the initial log likelihood
if  exitflag<=0 && LL<LL0
    %Try more iterations, bit only if the final likelihood is
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
    [parameters,LL,exitflag,output]=fmincon('egarch_likelihood',parameters,[],[],[],[],[],[],'egarch_nlcon',options,data_augmented,p,o,q,error_type,back_cast,T,1);
end


%If the optimization still hasn't converged, try other starting values
if startingflag==0 && exitflag<=0
    warning('Not Sucessful!  Trying alternative starting values.')
    %Keep track of the final estimates, if nothing converges, we will
    %return the best outcome and inform the user
    robust_parameters(1,:)=parameters';
    %Also keep the LL
    robust_LL(1)=LL;
    %And the output
    robust_output{1}=output;
    %Keep track of the iteration
    robust_iter=2;
    while exitflag<=0
        %This condition checks that we haven't converged
        %OR that the best objective is worse than best grid search
        %Sort the original grid search log likelihoods and parameters
        startingvals=ordered_parameters(robust_iter,:)';
        startingvals=[startingvals;nu;lambda];
        startingvals=egarch_transform(startingvals,p,o,q,error_type);

        LL0=egarch_likelihood(startingvals,data_augmented,p,o,q,error_type,back_cast,T,1);
        options.HessUpdate='bfgs';
        %Try the second set of starting values
        [parameters,LL,exitflag,output]=fmincon('egarch_likelihood',startingvals,[],[],[],[],[],[],'egarch_nlcon',options,data_augmented,p,o,q,error_type,back_cast,T,1);
        if  exitflag<=0 && LL<LL0
            %Again, if the LL improved, try more iterations
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
            [parameters,LL,exitflag,output]=fmincon('egarch_likelihood',parameters,[],[],[],[],[],[],'egarch_nlcon',options,data_augmented,p,o,q,error_type,back_cast,T,1);
        end
        %Save the parameter estimates, LL and output
        robust_parameters(robust_iter,:)=parameters';
        robust_LL(robust_iter)=LL;
        robust_output{robust_iter}=output;

        %Increment the robust_iter
        robust_iter=robust_iter+1;

        if robust_iter>size(ordered_parameters,1);
            %save the best LL and parameters and break
            warning('Convergance not achieved.  Use results with CAUTION.  You may need to provide starting values.');
            [LL,robust_iter]=min(robust_LL);
            parameters=robust_parameters(robust_iter,:)';
            output=robust_output{robust_iter};
            break
        end
    end
end

%Transform the parameters back to the usual space since the optimization is
%finished
parameters=egarch_itransform(parameters,p,o,q,error_type);

%Compute the log likelihood, the individal log likelihoods and the variances if needed
if nargout>1
    [LL, likelihoods, ht]=egarch_likelihood(parameters,data_augmented,p,o,q,error_type,back_cast,T);
    LL=-LL;
end

%Compute standard errors using RobustVCV if needed.
if nargout>3
    nw=0; %No newey west
    [VCVrobust,A,B,scores,hess]=robustvcv('egarch_likelihood',parameters,nw,data_augmented,p,o,q,error_type,back_cast,T);
    VCV=hess^(-1)/(T-m);
end

%Finally, set up the diagnostics.
diagnostics.EXITFLAG=exitflag;
diagnostics.ITERATIONS=output.iterations;
diagnostics.FUNCCOUNT=output.funcCount;
diagnostics.MESSAGE=output.message;
