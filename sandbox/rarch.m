function [parameters, ll, ht, VCV, scores, diagnostics] = rarch(data,p,q,type, method, startingvals,options)


% METHOD - [OPTIONAL] String, either 'Joint' or '2-step'.  Default is 2-step.
% TYPE - [OPTIONAL] String, either 'Scalar','Diagonal', or 'CP' (for common persistence). Default is
%        'Diagonal'


