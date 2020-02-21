function varargout = dirod(varargin)
% Date sorted folder listing
%
% USAGE:
%   dirod
%   dirod FILENAMEMASK
%   dirod('FILENAMEMASK')
%   LIST = dirod('FILENAMEMASK')
%
% INPUTS:
%   FILENAMEMASK - A file name mask, e.g. *, *.m or d*.m
%
% OUTPUTS:
%   LIST         - A directory structure ordered by date
%
% COMMENTS:
%   Similar to the windows command dir /od

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 6/7/2010

if isempty(varargin)
    arg = '*';
elseif length(varargin)==1
    arg = varargin{1};
else
    error('Too many input arguments.  0 or 1 input only');
end

dirData = dir(arg);
[~,index]=sort([dirData.datenum]);
dirData = dirData(index);

fileNameLen = length(dirData(1).name);
for i=2:length(dirData)
    fileNameLen = max(fileNameLen,length(dirData(i).name));
end

if nargout==0
    disp(' ');
    for i=1:length(dirData)
        disp([dirData(i).name repmat(' ',1,fileNameLen-length(dirData(i).name)+2) dirData(i).date]);
    end
    disp(' ');
else
    varargout{1} = dirData;
end
