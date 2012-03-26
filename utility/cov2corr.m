function [sigma,correl] = cov2corr(cov)

if ndims(cov)==2
    sigma = sqrt(diag(cov));
    correl = cov ./(sigma*sigma');
elseif ndims(cov)==3
    T = size(cov,3);
    K = size(cov,2);
    sigma = zeros(T,K);
    correl = zeros(K,K,T);
    for t=1:T
        sigma(t,:) = sqrt(diag(cov(:,:,t)))';
        correl(:,:,t) = cov(:,:,t)./(sigma(t,:)'*sigma(t,:));
    end
end