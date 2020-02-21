function [pac, pacstd, hfig] = spacf(data, lags, robust, graph)
% Computes sample partial autocorrelations and standard deviation using either heteroskedasticity
% robust standard errors or classic (homoskedastic) standard errors
%
% USAGE:
%  [PAC,PACSTD] = spacf(DATA,LAGS)
%  [PAC,PACSTD,HFIG] = spacf(DATA,LAGS,ROBUST)
%
% INPUTS:
%  DATA      - A T by 1 vector of data
%  LAGS      - The number of autocorrelations to compute
%  ROBUST    - [OPTIONAL] Logical variable (0 (non-robust) or 1 (robust)) to indicate whether
%                heteroskedasticity robust standard errors should be used. Default is to use robust
%                standard errors (ROBUST=1).
%  GRAPH     - [OPTIONAL] Logical variable (0 (no graph) or 1 (graph)) indicating whether the
%                function should produce a bar plot of the sample autocorrelations and confidence
%                intervals. Default is to produce a graphic (GRAPH=1).
%
% OUTPUTS:
%  PAC        - A LAGS by 1 vector of partial autocorrelations
%  PACSTD     - A LAGS by 1 vector of standard deviations
%  HFIG       - Figure handle to the bar plot of the autocorrelations
%
% COMMENTS:
%  Sample partial autocorrelations computed from autocorrelations that are computed using the
%  maximum number of observations for each lag.  For example, if DATA has 100 observations, the
%  first autocorrelation is computed using 99 data points, the second with 98 data points and so on.

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3.0.1    Date: 1/1/2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>4 || nargin<2
    error('2 to 4 inputs required.')
elseif nargin==2
    robust = true;
    graph = true;
elseif nargin==3
    graph = true;
end
T=length(data);
if T<=lags
    error('Length of data must be larger than LAGS')
end
if size(data,1)~=T,
    data=data';
end
if size(data,2)~=1
    error('DATA must be a column vector')
end
if ~isscalar(lags) || lags<0 || floor(lags)~=lags
    error('LAGS must be a positive integer')
end
if nargin==3
    if isempty(robust)
        robust=true;
    end
    if ~ismember(robust,[0,1])
        error('ROBUST must be either 0 or 1')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pac = zeros(lags,1);
pacstd = zeros(lags,1);
for i=1:lags
    [y, x] = newlagmatrix(data,i,1);
    t=length(y);
    rho = x\y;
    pac(i) = rho(i+1);
    if robust==1 && (graph || nargout > 1)
        e = y-x(:,1:i)*(x(:,1:i)\y);
        XpXi = ((x'*x)/t)^(-1);
        eX  = repmat(e,1,i+1).*x;
        XEX = (eX'*eX)/t;
        V = XpXi*XEX*XpXi/t;
        pacstd(i) = sqrt(V(i+1,i+1));
    else
        pacstd=ones(size(pac))/sqrt(length(data));
    end
end

if graph
    hfig=figure;
    set(hfig,'Position',[100 100 800 600])
    clf;
    h = bar(pac);
    hold on;
    h2 = plot((0:lags+1)',[2*[pacstd(1);pacstd;pacstd(lags)] -2*[pacstd(1);pacstd;pacstd(lags)]]);
    set(h,'FaceColor',[.5 .5 1]);
    axis tight;
    ax = axis;
    spread = .2*(ax(4)-ax(3));
    ax(1) = 0;
    ax(2) = lags+1;
    ax(3) = ax(3) - spread;
    ax(4) = ax(4) + spread;
    axis(ax);
    if robust
        t = title('Sample Partial Autocorrelations and Robust Standard Errors');
    else
        t = title('Sample Partial Autocorrelations and Non-robust Standard Errors');
    end
    set(t,'FontSize',14)
    set(h2(1),'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
    set(h2(2),'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
    hold off;
else
    hfig = [];
end
