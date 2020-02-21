function [autocorr, sigma2_y] = acf(phi,theta,N,sigma2_e)
% Computes the theoretical autocorrelations and long-run variance of an ARMA(p,q) process
%
% USAGE:
%  [AUTOCORR, SIGMA2_T] = acf(PHI,THETA,N)
%  [AUTOCORR, SIGMA2_T] = acf(PHI,THETA,N,SIGMA2_E)
% 
% INPUTS:
%  PHI         - Autoregressive parameters, in the order t-1,t-2,...
%  THETA       - Moving average parameters, in the order t-1,t-2,...
%  N           - Number of autocorrelations to be computed
%  SIGMA2_E    - [OPTIONAL] Variance of errors.  If omitted, sigma2_e=1
%
% OUTPUTS:
%  AUTOCORR    - N+1 by 1 vector of autocorrelation. To recover the autocovariance of an ARMA(P,Q),
%                  use AUTOCOV = AUTOCORR * SIGMA2_Y 
%  SIGMA2_Y    - Long-run variance, denoted gamma0 of ARMA process with innovation variance SIGMA2_E
%
% COMMENTS:
%  Note: The ARMA model is parameterized as follows:
%        y(t)=phi(1)y(t-1)+phi(2)y(t-2)+...+phi(p)y(t-p)+e(t)+theta(1)e(t-1)+theta(2)e(t-2)+...+theta(q)e(t-q)
%
%  To compute the autocorrelations for an ARMA that does not include all lags 1 to P, insert 0 for
%  any excluded lag.  For example, if the model was y(t) = phi(2)y(t-1) + e(t), PHI = [0 phi(2)] 

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 1/1/2007

if nargin<3 || nargin>4
    error('Requires 3 or 4 arguments.');
end

if length(N)>1 || any(N<0) || N~=floor(N)
    error('N must be a non-negative scalar integer.')
end

if nargin<4
    sigma2_e=1;
elseif length(sigma2_e)>1 || any(sigma2_e<=0)
    error('sigma2_e must be a positive scalar');
end

p=length(phi);
if min(size(phi))>1
    error('phi must be a vector');
end
if size(phi,1)~=p
    phi=phi';
end

q=length(theta);
if min(size(theta))>1
    error('theta must be a vector');
end
if size(theta,1)~=q
    theta=theta';
end

%The acf is only well defined for a stationary process, so we are going to
%check
if p>0
    [rho,stationary]=inverse_ar_roots(phi); %#ok<ASGLU>
    if ~stationary
        error('Autoregressive roots (phi) do not correspond to a stationary process');
    end
end

%save the orginal values
phi_original=phi;

%Ensure both have the same length
%Since we are working with a model of the form Phi(L)y(t)=Theta(L)e(t),
%we have to add a 1 and to reverse the sign of the autoregressive
%parameters
phi = [1 ; -phi; zeros(q-p,1) ];
theta = [1 ; theta; zeros(p-q,1)];
%m will be the dimension of our set of lienar equations, ti is the max(p,q)+1
m=max(p,q)+1;

%This is the left hand side of the set of linear equations
%These are all functions of the autoregressive parameters
%The autocovariance are found by solving a specific set of lienar equations
%of the form:
%phi_transformed gamma = theta_transformed delta
%where delta(i) is E[y_t e(t-i)], and the error variance is assumed to be 1
phi_transformed  = zeros(m,m);
T = toeplitz(1:m);
for i=1:m
    for j=1:m
        phi_transformed(i,T(i,j)) = phi_transformed(i,T(i,j))+phi(j);
    end
end

% On to the right hand side.
% First we need to calculate delta, which can be found by solving
% tril(toeplits(phi)) * delta = theta
delta=tril(toeplitz(phi))^(-1)*theta;

%The form of theta_transformed is
% [theta(0) theta(1)    theta(2)   ...         theta(q-1)  theta(q)]
% [theta(1) theta(2)    ...        theta(q-1)  theta(q)    0       ]
% [theta(2) theta(3)    ...        theta(q)    0           0       ]
% [...      ...         ...        ...         ...         ...     ]
% [theta(q) 0           0          ...         ...         0       ]
% This is accomplished by:
theta_transformed=flipud(tril(toeplitz(flipud(theta))));

%Finally, we can compute the first m autocovariances
autocov=phi_transformed^(-1)*(theta_transformed*delta);

%Finally, we need to return the correct number of autocovariances
%If fewer than m are required
if N+1<m
    autocov=autocov(1:N+1);
elseif N+1>m
    autocov=[autocov;zeros(N-m+1,1)];
    %Recursion to compute the autocovariances
    if p>0
        for i=m:N
            autocov(i+1)=phi_original'*autocov(i:-1:(i-p+1));
        end
    end
end
sigma2_y=sigma2_e*autocov(1);
autocorr=autocov./autocov(1);