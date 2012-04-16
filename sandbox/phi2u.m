function U = phi2u(phi)

m = length(phi);
k = ceil(sqrt(2*m));
U = eye(k);
c = cos(phi);
s = sin(phi);
count = 1;
for i=1:k
    for j=i+1:k
        R = eye(k);
        R(i,i) = c(count);
        R(i,j) = -s(count);
        R(j,i) = s(count);
        R(j,j) = c(count);
        U = U*R;
        count = count + 1;
    end
end
    