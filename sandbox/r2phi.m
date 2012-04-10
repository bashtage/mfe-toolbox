function phi = r2phi(R)

X = chol(R);
k = length(R);
S = zeros(k);

P = zeros(k);
S(1,:) = 1;
cumS = S;
for i=1:k-1
    P(i,i+1:k) = acos(X(i,i+1:k)./cumS(i,i+1:k));
    S(i+1,i+1:k) = sin(P(i,i+1:k));
    cumS(i+1,i+1:k) = prod(S(1:i+1,i+1:k));
end
phi = P(P>0);