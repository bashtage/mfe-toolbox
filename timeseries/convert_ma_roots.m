function parameters = convert_ma_roots(parameters, q)
% Checks an MA parameterization for invertability and inverts it if possible
%
%
%
% COMMENTS:
% It is not always possible to invert irregular MAs since the indices of the non-zero coefficients
% in the invertable MA may differ.  If this is the case PARAMETERS will be equal to the input.

MAparameters = zeros(1,max(q));
MAparameters(q) = parameters;
lambda = roots([1 MAparameters]);
lambda(abs(lambda)>1) = 1./lambda(abs(lambda)>1);
parameters = poly(lambda);
parameters(abs(parameters)<1e-10) = 0;

if length(q)<max(q)
    % Check for index change
    newq = find(abs(parameters)>1e-10);
    newq = newq(2:length(newq))-1;
    if length(q)~=length(newq) || any(q~=newq)
        warning('MFEToolbox:Invertability',['The irregular MA cannot be inverted without changing the lag indices. \nThe require lag indices are:' num2str(newq)])
        parameters = MAparameters(q);
    else
        parameters = parameters(q+1);
    end
else
    parameters = parameters(q+1);
end

if size(parameters,2)>size(parameters,1)
    parameters = parameters';
end