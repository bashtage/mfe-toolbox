function [rv,rvSS]=realized_threshold_multipower_variation(price,time,timeType,samplingType,samplingInterval,gamma,thresholdScale,thresholdType,subsamples)


c = thresholdScale;


logPrice =log(price);
% Filter prices and compute the RV
filteredLogPrice = realized_price_filter(logPrice,time,timeType,samplingType,samplingInterval);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filteredLogPrice  = log(price);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





returns = diff(filteredLogPrice);

L = 25;
V = inf*ones(size(returns));
K = -L:L;
K = 1/sqrt(2*pi)*exp(-(K/L).^2/2);
K(L+(-1:1)) = 0;
m = size(returns,1);
returns2 = returns.^2;

finished = false;
while ~finished 
    ind = returns2<(c^2.*V);
    Vold  = V;
    for i=1:m
        pl = i-L:i+L;
        valid = pl>1 & pl<=m;
        tempInd = ind(pl(valid));
        w = K(valid);
        V(i) = w*(returns2(pl(valid)).*tempInd)/(w*tempInd);
    end
    if all(V==Vold)
        finished = true;
    end
end

expectedValue = 1./(2*normcdf(-c)*sqrt(pi))*(2/c^2) * gamma(3/2).*gammainc(c^2/2,3/2,'upper') * c^2 *V;


rv = returns(ind)'*returns(ind) + sum(expectedValue(~ind));