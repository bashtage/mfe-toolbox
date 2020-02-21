function [varargout]  = tsresidualplot(y,errors,dates)
% Produces a plot for visualizing time series data and residuals from a time series model
%
% USAGE:
%   tsresidualplot(Y,ERRORS)
%   [HAXIS,HFIG] = tsresidualplot(Y,ERRORS,DATES)
%
% INPUTS:
%   Y      - A T by 1 vector of data
%   ERRORS - A T by 1 vector of residuals, usually produced by ARMAXFILTER
%   DATES  - [OPTIONAL] A T by 1 vector of MATLAB dates (i.e. should be 733043 rather than '1-1-2007').
%              If provided, the data and residuals will be plotted against the date rather than the
%              observation index
%
% OUTPUTS:
%   HAXIS  - A 2 by 1 vector axis handles to the top subplots
%   HFIG   - A scalar containing the figure handle
%
% COMMENTS:
%   HAXIS can be used to change the format of the dates on the x-axis when MATLAB dates are provides
%   by calling 
%
%        datetick(HAXIS(j),'x',DATEFORMAT,'keeplimits')
%
%   where j is 1 (top) or 2 (bottom subplot) and DATEFORMAT is a numeric value between 28.  See doc
%   datetick for more details.  For example, 
%
%        datetick(HAXIS(1),'x',25,'keeplimits')
%
%   will change the top subplot's x-axis labels to the form yy/mm/dd.
%
% EXAMPLES:
%   Estimate a model and produce a plot of fitted and residuals
%       [parameters, LL, errors] = armaxfilter(y, 1, 1, 1);
%       tsresidualplot(y, errors)
%   Estimate a model and produce a plot of fitted and residuals with dates
%       [parameters, LL, errors] = armaxfilter(y, 1, 1, 1);
%       dates = datenum('01Jan2007') + (1:length(y));
%       tsresidualplot(y, errors, dates)
%
% See also ARMAXFILTER, DATETICK

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 1/1/2007


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>3 || nargin<2
    error('2 or 3 inputs required.')
end
T=length(y);
if size(y,2)>size(y,1)
    y=y';
end
if size(errors,2)>size(errors,1)
    errors=errors';
end

if size(y,2)~=1 || ndims(y)~=2
    error('Y must be a column vector of data')
end

if size(errors,2)~=1 || ndims(errors)~=2
    error('Y must be a column vector of data')
end
if ~all(size(y)==size(errors))
    error('Y and ERRORS must have the same dimensions')
end

if nargin<3 || isempty(dates)
    dates=(1:T)';
    dateflag=0;
else
    if size(dates,2)>size(dates,1)
        dates=dates';
    end
    if size(dates,2)~=1 || ndims(dates)~=2
        error('DATES must be a column vector with the same dimensions as Y')
    end
    if ndims(dates)~=2
        error('DATES must be a column vector with the same dimensions as Y')
    end
    if ~isnumeric(dates)
        error('DATES must contain numeric MATLAB dates, not strings.')
    end
    dateflag=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


yhat = y-errors;
% Fit and Actual subplot
fig = figure('Position',[100 100 800 600]);
s1=subplot(2,1,1);
h1=plot(dates,[y yhat]);
axis tight;
ax=axis;
spread = ax(4)-ax(3);
ax(3)=ax(3)-.1*spread;
ax(4)=ax(4)+.1*spread;
axis(ax);
legend('Data','Fit');
if dateflag
    datetick('x','keeplimits')
end
t1=title('Data and Fit');
set(h1,'LineWidth',2)
set(s1,'LineWidth',2,'FontSize',12,'FontWeight','Bold')
set(t1,'FontSize',12,'FontWeight','Bold')
% Residual Subplot
s2=subplot(2,1,2);
h2=plot(dates,errors);
axis tight;
ax=axis;
spread = ax(4)-ax(3);
ax(3)=ax(3)-.1*spread;
ax(4)=ax(4)+.1*spread;
axis(ax);
legend('Residual');
if dateflag
    datetick('x','keeplimits')
end
t2=title('Residual');
set(h2,'LineWidth',2)
set(s2,'LineWidth',2,'FontSize',12,'FontWeight','Bold')
set(t2,'FontSize',12,'FontWeight','Bold')

switch nargout
    case 0
    case 1
        varargout{1}=[s1;s2];
    case 2
        varargout{1}=[s1;s2];
        varargout{2} = fig;
    otherwise
        error('0, 1 or 2 output arugments supported.')
end
