function [ll,lls,Rt] = dcc_likelihood(parameters,data,dataAsym,m,l,n,R,N,backCast,backCastAsym,stage,transformIntercept,univariate)


[k,~,T] = size(data);
offset = 0;
H = zeros(T,k);
if stage==1
    for i=1:k
        u = univariate{i};
        count = u.p+u.o+u.q+1;
        parameters = parameters(offset + (1:count));
        offset = offset+count;
        ht = tarch_core(u.fdata,u.fIdata,parameters,u.back_cast,u.p,u.o,u.q,u.m,u.T,u.tarch_type);
        H(:,i) = ht(m+1:T);
    end
    stdData = zeros(k,k,T);
    stdDataAsym = zeros(k,k,T);
    for t=1:T
        h = sqrt(H(t,:));
        stdData(:,:,t) = data./(h*h');
        stdDataAsym(:,:,t) = dataAsym./(h*h');
    end
end
if stage<=2
    z = parameters(offset + (1:k*(k-1)/2));
    offset = k*(k-1)/2;
    if transformIntercept
        R = z2r(z);
    else
        R = corr_ivech(z);
    end
end
a = parameters(offset + (1:m));
g = parameters(offset + (m+1:l));
b = parameters(offset + (m+l+1:m+l+n));

if stage>=2
    stdData = data;
    stdDataAsym = dataAsym;
    logdetH = zeros(T,1);
else
    logdetH = sum(log(H),2);
end

if stage==3
    intercept = R*(1-sum(a)-sum(b))-N*sum(g);
else 
    intercept = R*(1-sum(a)-sum(b)-0.5*sum(g));
end
% Check eigenvalues?

I = eye(k);
Qt = zeros(k,k,T);
Rt = zeros(k,k,T);
lls = zeros(T,1);
for t=1:T
    Qt(:,:,t) = intercept;
    for i = 1:m
        if (t-i)>0
            Qt(:,:,t) = Qt(:,:,t) + a(i)*stdData(:,:,t-i);
        else
            Qt(:,:,t) = Qt(:,:,t) + a(i)*backCast;
        end
    end
    for i = 1:l
        if (t-i)>0
            Qt(:,:,t) = Qt(:,:,t) + g(i)*stdDataAsym(:,:,t-i);
        else
            Qt(:,:,t) = Qt(:,:,t) + g(i)*backCastAsym;
        end
    end
    for i = 1:n
        if (t-i)>0
            Qt(:,:,t) = Qt(:,:,t) + b(i)*Qt(:,:,t-i);
        else
            Qt(:,:,t) = Qt(:,:,t) + b(i)*backCast;
        end
    end
    n = sqrt(diag(Q(:,:,t)));
    Rt(:,:,t) = Qt(:,:,t)./ (n*n');
    lls(t) = 0.5*(likconst + logdetH(t) + log(det(Rt(:,:,t))) + sum(diag((Rt\I)*stdData)));
end
ll = sum(lls);

if isnan(ll) || ~isreal(ll) || isinf(ll)
    ll = 1e7;
end