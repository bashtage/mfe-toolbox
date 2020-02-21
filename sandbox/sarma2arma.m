function [parameters, arP, maQ] = sarma2arma(parameters, p, q, seasonal)
% Parameter order: c, p, Sp(1),Sp(2),...,Sp(K),q, Sq(1),Sq(2),...,Sq(K)

Sp = seasonal(:,1);
Sq = seasonal(:,3);
S = seasonal(:,4);

arParamters = parameters(1:length(p)+sum(Sp));
[arParameters,arP]=convertSARtoAR(arParamters, p, Sp, S);
maParamters = parameters(length(p)+sum(Sp)+1:length(p)+sum(Sp)+length(q)+sum(Sq));
[maParameters,maQ]=convertSARtoAR(maParamters, q, Sq, S);
parameters = [arParameters maParameters]';

function [parameters, lags]=convertSARtoAR(parameters,p,Sp,S)
offset = 0;
if isempty(p) 
    poly = 1;
else
    lagPoly = zeros(1,max(p)+1);
    lagPoly(1) = 1;
    lagPoly(p+1) = -parameters(offset+1:offset+length(p));
    offset = offset + length(p);
    poly = lagPoly;
end
for i=1:length(S)
    p = Sp(i);
    sLag = S(i);
    if p>0
        lagPoly = zeros(1,1+sLag*p);
        lagPoly(1) = 1;
        lagPoly(sLag+1:sLag:(sLag*p+1)) = -parameters(offset+1:offset+p);
        offset = offset + p;
        poly = conv(poly, lagPoly);
    end
end
poly = poly(2:length(poly));
lags = find(poly~=0);
parameters = -poly(lags);



