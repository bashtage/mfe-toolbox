% This m-file produces the scales used in realized_range.  It can be run with more simulations to
% increase the precision of the estimated scales used.  If you do this you will have to copy and
% paste the results into the function since it uses a lookup table and interpolation.
 
% Clear 
clear all
% Number of simulations
BB = 1000000;
 
% Number of prices to use in each interval.  All integer factors of 23400.  The asymptotic value is
% 4*log(2)
 
m1 = factor(23400);
m1 = [1 m1];
ms = m1';
for i=2:length(m1);
    temp = nchoosek(1:length(m1),i);
    ms = [ms;    unique(prod(m1(temp),2))]; %#ok<AGROW>
end
ms = unique(ms);
maxM = max(ms);
 
% Turn the number of prices per window in to a step size
gap = maxM./ms;
ms = ms+1;
% Initialize a place to hold results
M = length(ms);
MC = zeros(BB,M);
% BB simultions
tic
for j = 1:BB
    % Resample a single BM 1, 2, ..., 23400 times.
    x = [0;cumsum(randn(maxM,1)/sqrt(maxM))];
    for i=1:M
        x2 = x(1:gap(i):(maxM+1));
        MC(j,i)  = (max(x2)-min(x2))^2;
    end
    % Display the count and the expected time remaining.
    if mod(j,10000)==0
        t=toc;
        str = [num2str(j) ' iterations complete.  The elapsed time is ' num2str(t) '. The expected time remaining is ' num2str(BB*t/j) ];
        disp(str)
    end
end
 
 
% Some smoothing using a concave regression
y = mean(MC);
y(1) = 1;
orig_mean_MC=y;
n = length(y);
x = eye(n);
A1 = [eye(n-1) zeros(n-1,1)] + [zeros(n-1,1) -eye(n-1)];
A2 = [eye(n-2) zeros(n-2,2)] + [zeros(n-2,1) -2*eye(n-2) zeros(n-2,1)] + [zeros(n-2,2) eye(n-2)];
b = zeros(2*n-3,1);
concave_mean_MC=lsqlin(x,y,[A1;A2],b);
% Save the results.  
%save realized_range_simulation_results orig_mean_MC concave_mean_MC ms

