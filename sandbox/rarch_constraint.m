function [c,ceq] = rarch_constraint(parameters,data,p,q,C,backCast,type,isJoint)

ceq = [];
k = size(data,1);
[~,A,B] = rarch_parameter_transform(parameters,p,q,k,C,type,isJoint);
constraint = diag(sum(A.^2,3)+sum(B.^2,3) - .99998);
switch type
    case 1
        c = constraint(1);
    case 2
        c = constraint;
    case 3
        c = constraint;
end