function [pautocorr] = pacf(phi,theta,N)
% Computes the theoretical partial autocorrelations an ARMA(p,q) process
%
% USAGE:
%  [PAUTOCORR] = pacf(PHI,THETA,N)
% 
% INPUTS:
%  PHI       - Autoregressive parameters, in the order t-1,t-2,...
%  THETA     - Moving average parameters, in the order t-1,t-2,...
%  N         - Number of autocorrelations to be computed
%
% OUTPUTS:
%  PAUTOCORR - N+1 by 1 vector of partial autocorrelations. 
%
% COMMENTS:
%  Note: The ARMA model is parameterized:
%        y(t)=   phi(1)y(t-1) +   phi(2)y(t-2) + ... +   phi(p)y(t-p)
%             +theta(1)e(t-1) + theta(2)e(t-2) + ... + theta(q)e(t-q) + e(t)
%
%  To compute the autocorrelations for an ARMA that does not include all lags 1 to P, insert 0 for
%  any excluded lag.  For example, if the model was y(t) = phi(2)y(t-1) + e(t), PHI = [0 phi(2)] 

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 1/1/2007


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin~=3
    error('3 inputs required.')
end
if ~isscalar(N) && N>0 && floor(N)==N
    error('N must be a positive integer')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ac = acf(phi,theta,N+1);
ac = ac(2:N+1);

N=N;
pac = zeros(N,1);
% Begin the recursion by computing the first partial autocorrelation, which
% is simply the first autocorrelation
pac(1) = ac(1);
% This code uses a partitioned inverse to compute the PACF
if N>=2
    % The the second, here we will use regression formulas
    XpX = toeplitz([1 ac(1)]);
    XpXinv = XpX^(-1);
    Xpy = ac(1:2);
    temp = XpXinv * Xpy;
    pac(2) = temp(2);
    % Now to compute the remaining using the partitioned inverse
    for i=3:N
        Ainv = XpXinv;
        B = ac(i-1:-1:1);
        C = B';
        D = 1;
        SDinv = Ainv + Ainv*B*(D-C*Ainv*B)^(-1)*C*Ainv;

        XpXinv = [SDinv -SDinv*B*D^(-1);
            -D^(-1)*C*SDinv D^(-1)+D^(-1)*C*SDinv*B*D^(-1)];
        XpXinv = (XpXinv+XpXinv')/2;
        Xpy = ac(1:i);
        temp = XpXinv * Xpy;
        pac(i) = temp(i);
    end
end

pautocorr = [ 1;pac];
pautocorr(abs(pautocorr)<100*eps)=0;