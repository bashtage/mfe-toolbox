function z = r2z(R)

k = length(R);
C = zeros(k);
C2 = chol(R)';

for i=2:k
    rem = 1;
    for j=i-1:-1:1
        C(i,j) = C2(i,j)/sqrt(rem);
        rem = rem - C2(i,j)^2;
    end
end

C =C';
z = C(~tril(true(k)));
z = log((z+1)./(1-z));