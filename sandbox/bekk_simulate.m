function [data, Ht] = bekk_simulate(T,k,parameters,p,o,q,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isscalar(T)
    e = randn(2*T,k);
else
    e = T;
    if size(e,2)~=k
        error('T must have K columns when providing simulated random numbers.')
    end
    T = size(e,1);
    e = [e(ceil(rand(T,1)*T),:);e];
end

if strcmpi(type,'Scalar')
    type = 1;
elseif strcmpi(type,'Diagonal')
    type = 2;
elseif strcmpi(type,'Full')
    type = 3;
else
    error('TYPE must be ''Scalar'', ''Diagonal'' or ''Full''.')
end

if p<=0 || floor(p)~=p
    error('P must be a positive scalar integer.')
end
if o<0 || floor(o)~=o
    error('O must be a positive scalar integer.')
end
if q<0 || floor(q)~=q
    error('Q must be a positive scalar integer.')
end

k2 = k*(k+1)/2;
switch type
    case 1
        count = p+o+q;
    case 2
        count = (p+o+q)*k;
    case 3
        count = (p+o+q)*k*k;
end
count = count + k2;
if length(parameters)~=count
    error('PARAMETERS does not have the expected number of elements.')
end
[C,A,G,B] = bekk_parameter_transform(parameters,p,o,q,k,type);
m = zeros(k*k);
for i=1:p
    m = m + kron(A(:,:,i),A(:,:,i));
end
for i=1:o
    m = m + 0.5*kron(G(:,:,i),G(:,:,i));
end
for i=1:q
    m = m + kron(B(:,:,i),B(:,:,i));
end

if max(eigs(m))>1
    backCast = C/.001;
    warning('MFE:nonstationary','The parameters do not correspond to the stationary region.')
else
    backCast = inv(eye(k*k)-m)*C(:);
    backCast = reshape(backCast,k,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ht = repmat(eye(k),[1 1 2*T]);
eta = zeros(size(e));
for i=1:(2*T)
    Ht(:,:,i) = C;
    for j=1:p
        if (i-j)<=0
            Ht(:,:,i) = Ht(:,:,i) + A(:,:,j)'*backCast*A(:,:,j);
        else
            Ht(:,:,i) = Ht(:,:,i) + A(:,:,j)'*(e(i-j,:)'*e(i-j,:))*A(:,:,j);
        end
    end
    for j=1:o
        if (i-j)<=0
            Ht(:,:,i) = Ht(:,:,i) + G(:,:,j)'*backCast*G(:,:,j);
        else
            Ht(:,:,i) = Ht(:,:,i) + G(:,:,j)'*(eta(i-j,:)'*eta(i-j,:))*G(:,:,j);
        end
    end    
    for j=1:q
        if (i-j)<=0
            Ht(:,:,i) = Ht(:,:,i) + B(:,:,j)'*backCast*B(:,:,j);
        else
            Ht(:,:,i) = Ht(:,:,i) + B(:,:,j)'*Ht(:,:,i-j)*B(:,:,j);
        end
    end
    Ht12 = Ht(:,:,i)^(0.5);
    e(i,:) = e(i,:)*Ht12;
    eta(i,:) = e(i,:).*(e(i,:)<0);
end

data = e(T+1:2*T,:);
Ht = Ht(:,:,T+1:2*T);