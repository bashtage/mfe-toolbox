function fr = weights_to_frequency_response(a,w,n)
% Computes the frequency response on [0,0.5] for a transformation of the form
%
%      A'y = W'x
%
% where A is a K by 1 vector and W is an M by 1 vector.  
%
% USAGE:
%   [fr] = weights_to_frequency_response(A,W,N)
%
% INPUTS:
%   A  - K by 1 vector of output weights
%   W  - M by 1 vector of input weights
%   N  - Number of points to use when computing the frequency response (e.g. linspace(0,0.5,N))
%
% OUTPUTS:
%   FR - N by 1 vector containing the gain
%
% COMMENTS:
%
% EXAMPLES:
%   Centered MA(3)
%       fr =  weights_to_frequency_response([0 1 0]',[1/3 1/3 1/3])

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 4    Date: 11/1/2009

if nargin==2
    n = 100;
end
    
f = linspace(0,0.5,n);
if size(w,2)>size(w,1)
    w=w';
end
if size(a,2)>size(a,1)
    a=a';
end
m = length(w);
k = length(a);
fr = zeros(size(f));
for j=1:length(f)
    fr(j) = (w'*exp(-1i*(1:m)'*2*pi*f(j)))/(a'*exp(-1i*(1:k)'*2*pi*f(j)));
end
fr = abs(fr);   
