function [c,ceq]=egarch_nlcon(parameters, data, p, o, q, error_type, back_cast, T, estim_flag)
%  EGARCH(P,O,Q) nonlinear parameter restriction.  Used in estimation of
%  EGARCH models.
%
%  USAGE:
%    [C,CEQ]=egarch_nlcon(PARAMETERS,DATA,P,O,Q,ERROR_TYPE,BACK_CAST,T,ESTIM_FLAG)
% 
%  INPUTS:
%    See EGARCH_LIKELIHOOD
% 
%  OUTPUTS:
%    C    - Vector of nonlinear inequality constraints.  Based on the roots
%           of a polynomial in beta
%    CEQ  - Empty matrix
% 
%  COMMENTS:
% 
%  See also EGARCH

beta=parameters(p+o+2:p+o+q+1);
%zn=fliplr(-beta');
%Find any zeros
%if any(abs(zn)<1e-08);
%    zn(abs(zn)<1e-08)=sign(zn(abs(zn)<1e-08))*1e-08;
%end
%r=roots([zn 1]);
%c=1./abs(r)-.9998;
c=abs(roots([1;-beta]))-.99998;
ceq=[];