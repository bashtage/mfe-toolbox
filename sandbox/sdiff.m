function [y,lags,scales] = sdiff(x, d, S)

poly = 1;
for i=1:length(d)
    for j=1:d(i)
        poly = conv(poly,[1 zeros(1,S(i)-1) -1]);
    end
end
lags = find(poly~=0);
n = length(lags);
lags = lags(2:n);
scales = poly(lags);
lags = lags - 1;
y = [];

if ~isempty(x)
    ml = max(lags);
    T = length(x);
    if (ml+1)<=T
        y = x(ml+1:T);
        for j=2:length(lags)
            y = y + scales(j) * x(ml+1-lags(j):T-lags(j));
        end
    end
end