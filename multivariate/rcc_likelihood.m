function [ll,lls,Rt] = rcc_likelihood(parameters,data,m,n,R,backCast,stage,type,composite,isJoint,isInference,rScale,univariate)

% Parse parameters
% Type, 1: Scalar, 2: CP, 3:Diagonal
% Stage: 1: Joint, 2: Intercepts and dynamcs, 3: Dynamics
% IsInference: True: no transform of intercept, False: transform intercept
% IsJoint: True: Parameters has vol, intercept and dynamics

[k,~,T] = size(data);
offset = 0;
% Parse Parameters
if stage==1 || isJoint
    count = 0;
    for t=1:k
        u = univariate{t};
        count = count + u.p+u.o+u.q+1;
    end
    garchParameters = parameters(1:count);
    offset = offset + count;
    computeVol = true;
else
    computeVol = false;
end
if stage<=2 || isJoint
    count = k*(k-1)/2;
    R = parameters(offset + (1:count));
    offset = offset + count;
    if isInference
        R = corr_ivech(R);
    else
        R = z2r(R);
    end
end
% Use rarch_parameter_transform with isJoint=false
[R,A,B] = rarch_parameter_transform(parameters(offset + (1:length(parameters)-offset)),m,n,k,R,type,false,false);
% Fix B in case of CP
B(B<0)= 0;

% Compute volatilities
H = ones(T,k);
if computeVol
    H = dcc_reconstruct_variance(garchParameters,univariate);
    stdData = zeros(k,k,T);
    for t=1:T
        h = sqrt(H(t,:));
        stdData(:,:,t) = data(:,:,t)./(h'*h);
    end
else
    stdData = data;
end


Qt = zeros(k,k,T);
e = zeros(k,k,T);
lls = zeros(T,1);

R = R .* sqrt(rScale*rScale');
R12 = R^(0.5);
Rm12 = R^(-0.5);
intercept = eye(k) - sum(A.^2,3) - sum(B.^2,3);
dint = diag(intercept);
dint(dint<.000001)=.000001;
intercept = diag(dint);

% Indices or constant, as needed
if composite == 0
    likconst = k*log(2*pi);
elseif composite == 1
    indices = [(1:k-1)' (2:k)'];
elseif composite == 2
    [i,j] = meshgrid(1:k);
    indices = [i(~triu(true(k))) j(~triu(true(k)))];
end

Rt = zeros(k,k,T);
for t=1:T
    e(:,:,t) = Rm12 * stdData(:,:,t) * Rm12;
    Qt(:,:,t+1) = intercept;
    for j=1:m
        if (t-j)>0
            Qt(:,:,t) = Qt(:,:,t) + A(:,:,j)*e(:,:,t-j)*A(:,:,j);            
        else
            Qt(:,:,t) = Qt(:,:,t) + A(:,:,j)*backCast*A(:,:,j);
        end
    end
    for j=1:n
        if (t-j)>0
            Qt(:,:,t) = Qt(:,:,t) + B(:,:,j)*Qt(:,:,t-j)*B(:,:,j);
        else
            Qt(:,:,t) = Qt(:,:,t) + B(:,:,j)*backCast*B(:,:,j);
        end
    end
    R = R12*Qt(:,:,t)*R12;
    rr = sqrt(diag(R)*diag(R)');
    R = R ./ rr;
    Rt(:,:,t) = R;
    h = sqrt(H(t,:));
    hh = h'*h;
    V = R.*hh;
    V = (V + V')/2;
    if composite == 0
        lls(t) = 0.5*(likconst + log(det(V)) + sum(diag(V^(-1)*data(:,:,t))));
    elseif composite
        lls(t) = composite_likelihood(V,data(:,:,t),indices);
    end
end
ll = sum(lls);