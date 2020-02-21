clear all
clc
close all
parameters = [.3 .05 .2 .4 .7 .55];
p = [0 1;0 1];
q = eye(2);
m = [1 390];
T = 1000;
[data,htOrig] = heavy_simulate(T,2,parameters,p,q,m);
parametersOrig=parameters;
data2 = data;
data2(:,1) = data2(:,1).^2;
backCast = mean(data2);
lb = min(data2)'/10000;
ub = max(data2)'*100000;
%[ll,lls,h] = heavy_likelihood(parameters,data2,p,q,backCast,lb,ub);

options = optimset('fminunc');
options.Display = 'iter';

sv = parameters;
sv(1:2) = log(sv(1:2));
%out = fminunc(@heavy_likelihood,sv,options,data2,p,q,backCast,lb,ub);

[parameters, ll, ht, VCV, scores] = heavy(data,p,q,'None');
plot([ht(:,1),htOrig(:,1)])


[parameters, ll, ht, VCV, scores] = heavy(data,p,q,'None',parametersOrig');
