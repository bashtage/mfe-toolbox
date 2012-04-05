function R = ZtoR(z)

m = length(z);
k = ceil(sqrt(2*m));
if k*(k-1)/2~=m
    error('Incorrect number of elements in z')
end
C = zeros(k);
z = (exp(z)-1)./(1+exp(z));
count = 1;
for i=2:k
    for j=1:i-1;
        C(i,j) = z(count);
        count = count + 1;
    end
end
C(1,1) = 1;
for i=2:k
    rem = 1;
    for j=i-1:-1:1
        C(i,j) = C(i,j)*sqrt(rem);
        rem = rem - C(i,j)^2;
    end
    C(i,i) = sqrt(rem);
end

R = C*C';
r = sqrt(diag(R));
R = R ./ (r*r')