% This m-file computes the scales required for realized_quantile_variance for samples with
% m={5 6 10 13 15 18 20 25 26 30 36 39 50 60 65 72 75 100 144} returns.  This substantially
% increases the speed fo realized quantile variance.  If using realized_quantile_variance with some
% number of returns per interval not listed, you can add it to the list and re-run this program.  It
% is also necessary to update one line in realized_quantile_variance.
% Clear everything
clear all

MS=[2 3 4 5 6 8 9 10 12 13 15 18 20 24 25 26 30 36 39 40 50 60 65 72 75 100 144];
% Initialize the place holders for the output
symmetricSimulationSamplerperbin = cell(2,1);
symmetricSimulationQuantile = cell(2,1);
asymmetricSimulationSamplerperbin = cell(2,1);
asymmetricSimulationQuantile = cell(2,1);
symmetricExpectedCovariance = cell(2,1);
symmetricExpectedQuantiles = cell(2,1);
asymmetricExpectedCovariance = cell(2,1);
asymmetricExpectedQuantiles = cell(2,1);
% Count is used since M is irregular
count = 1;
sumM = 0;
% Set the number of simulations
simulations = 100000000;
% Initialize the timer

options = optimset('lsqlin');
options.LargeScale = 'off';
% Symmetric
tic;
for M=[2 3 4 5 6 8 9 10 12 13 15 18 20 24 25 26 30 36 39 40 50 60 65 72 75 100 144];
    disp('Symmetric');
    disp('M=')
    disp(M)
    disp(' ');
    % Core calls realized_quantile_variance_scale
    symmetricSimulationSamplerperbin{count} = M;
    symmetricSimulationQuantile{count} = (1:M)/M;
    
    [symmetricExpectedQuantiles{count},temp,symmetricExpectedCovariance{count}]=realized_quantile_variance_scale(M,symmetricSimulationQuantile{count},simulations,true);
    sumM = sumM + M;
    disp(['Time elapsed is ' num2str(toc) '.']);
    
    
    if M>=10
        % Convex regression for some smoothing when M>=10
        y = symmetricExpectedQuantiles{count}';
        X=eye(length(y));
        n = length(y);
        A=-[eye(n-2) zeros(n-2,2)] + [zeros(n-2,1) 2*eye(n-2) zeros(n-2,1)] -[zeros(n-2,2) eye(n-2)];
        b = zeros(n-2,1);
        symmetricExpectedQuantiles{count} = lsqlin(X,y,A,b,[],[],[],[],[],options)';
        disp(max(abs(symmetricExpectedQuantiles{count}-y')));
    end
    count = count + 1;
    % Save each iteration
    save realized_quantile_scales symmetric* 
end
% Save
save realized_quantile_scales symmetric*




count = 1;
sumM = 0;
% Set the samples per interval that will be used
for M=[4 5 6 10 13 15 18 20 25 26 30 36 39 50 60 65 72 75 100 144];
    disp('Asymmetric ');
    disp('M=')
    disp(M)
    disp(' ');
    % Core calls realized_quantile_variance_scale
    asymmetricSimulationSamplerperbin{count} = M;
    asymmetricSimulationQuantile{count} = (ceil(M/2)+1:M)/M;
    [asymmetricExpectedQuantiles{count},temp,asymmetricExpectedCovariance{count}]=realized_quantile_variance_scale(M,asymmetricSimulationQuantile{count},simulations );
    sumM = sumM + M;
    disp(['Time elapsed is ' num2str(toc) '. '])
    
    
    if M>=10
        % Convex regression for some smoothing when M>=10
        y = asymmetricExpectedQuantiles{count}';
        X=eye(length(y));
        n = length(y);
        A=-[eye(n-2) zeros(n-2,2)] + [zeros(n-2,1) 2*eye(n-2) zeros(n-2,1)] -[zeros(n-2,2) eye(n-2)];
        b = zeros(n-2,1);
        asymmetricExpectedQuantiles{count} = lsqlin(X,y,A,b,[],[],[],[],[],options)';
        disp(max(abs(asymmetricExpectedQuantiles{count}-y')));
    end
    count = count + 1;
    % Save each iteration
    save realized_quantile_scales symmetric* asymmetric*
end
% Save
save realized_quantile_scales symmetric* asymmetric*
