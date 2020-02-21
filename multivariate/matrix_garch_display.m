function [C,A,G,B] = matrix_garch_display(parameters,p,o,q)
% Displays the MATRIX garch estimates in a friendly format
%
% USAGE:
%  [C,A,G,B] = matrix_garch(PARAMETERS,P,O,Q)
%
% INPUTS:
%   PARAMETERS   - Output from callig matrix_garch
%   P            - Positive, scalar integer representing the number of lags of the innovation process
%   O            - Non-negative scalar integer representing the number of asymmetric lags to include
%   Q            - Non-negative scalar integer representing the number of lags of conditional covariance
%
% OUTPUTS:
%   C           - K by K matrixd containing the constant
%   A           - A [K K P] dimension matrix containing the parameters of the symmetric terms
%   G           - A [K K O] dimension matrix containing the parameters of the asymmetric terms
%   B           - A [K K B] dimension matrix containing the parameters of the lagged conditional variance
%
% COMMENTS:
%
% EXAMPLES:
%     Display the parameters of an asymmetric MATRIX GARCH(1,1,1) process
%       [simulatedData, Ht, pseudoRC] = scalar_vt_vech_simulate(1000, [.02 .04 .95], .01*eye(2), 1, 1, 1, 72);
%       parameters = matrix_garch(simulatedData,[],1,1,1)
%       [C,A,G,B] = matrix_garch_display(parameters,1,1,1)

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 3/10/2011


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(parameters,2)>size(parameters,1)
    parameters = parameters';
end
if any(floor([p o q])~=[p o q]) || any ([p o q]<0)
    error('P, O, and Q must all be non-negative integers.')
end
% Dimension 
paraLen=length(parameters);
k2 = paraLen/(1+p+o+q);
if ~(floor(k2)==k2)
    error('PARAMETERS does not have the expected number of elements: K*(K+1)/2*(1+P+O+Q)')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


k = floor(sqrt(2*k2));


parameterMatrices= zeros(k,k,1+p+o+q);
index = 0;
for i=1:(1+p+o+1)
    temp = vec2chol(parameters(index+1:index+k2));
    parameterMatrices(:,:,i) = temp*temp';
    index=index+k2;
end

C = parameterMatrices(:,:,1);
A = parameterMatrices(:,:,2:p+1);
G = parameterMatrices(:,:,2+p:1+p+o);
B = parameterMatrices(:,:,2+p+o:1+p+o+q);

k1 = ceil(k/2);
dispStr = [repmat('      ',k1-1,1); 'H(t) =';repmat('      ',k-k1,1)];
dispStr = [dispStr repmat(' ',k,1) repmat('[',k,1) num2str(C) repmat(']',k,1)];
plus = [repmat(' ',k,1)  [repmat(' ',k1-1,1); '+';repmat(' ',k-k1,1)]];
for i=1:p
    dispStr = [dispStr plus repmat(' ',k,1) repmat('[',k,1) num2str(A(:,:,i)) repmat(']',k,1)]; %#ok<*AGROW>
    dispStr = [dispStr  [repmat(' ',k1-1,15); '.*r(t-' num2str(i) ')r(t-' num2str(i) ')''';repmat(' ',k-k1,15)]];
end
for i=1:o
    dispStr = [dispStr plus repmat(' ',k,1) repmat('[',k,1) num2str(G(:,:,i)) repmat(']',k,1)];
    dispStr = [dispStr  [repmat(' ',k1-1,15); '.*n(t-' num2str(i) ')n(t-' num2str(i) ')''';repmat(' ',k-k1,15)]];
end   
for i=1:q
    dispStr = [dispStr plus repmat(' ',k,1) repmat('[',k,1) num2str(B(:,:,i)) repmat(']',k,1)];
    dispStr = [dispStr  [repmat(' ',k1-1,8); '.*H(t-' num2str(i) ')';repmat(' ',k-k1,8)]];
end    

disp(dispStr)
if o>0
    disp('Note: r(t-i)r(t-i)'' is the RC, and n(t-i)n(t-i)'' is the asymmetric RC term, if the model used realized measures.')
else
    disp('Note: r(t-i)r(t-i)'' is the RC if the model used realized measures.')
end

