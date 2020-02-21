function addToPath(arg)

% Adds the MFE Toolbox to the current path and optionally saves the path
silent = false;
if nargin>0 && strcmpi(arg,'-silent')
    silent = true;
end

dirList = {'bootstrap',...
    'crosssection',...
    'distributions',...
    'GUI',...
    'multivariate',...
    'tests',...
    'timeseries',...
    'univariate',...
    'utility',...
    'realized',...
    };
optionalDirList = {'duplication'};
mexDir = 'dlls';

for i=1:length(dirList)
    addpath(fullfile(pwd,dirList{i}));
end


user_entry = '';
while ~ismember({user_entry},{'Y','N'}) && ~silent
    user_entry = input('Add MATLAB work-a-like functions? [Y/N(default)] ','s');
    user_entry = upper(user_entry);
    if isempty(user_entry)
        user_entry = 'N';
    end
end

if strcmpi(user_entry,'Y') 
    for i=1:length(optionalDirList)
        addpath(fullfile(pwd,optionalDirList{i}));
    end
end

user_entry = '';
while ~ismember({user_entry},{'Y','N'}) && ~silent
    user_entry = input('Make path permanent? [Y/N (default)] ','s');
    user_entry = upper(user_entry);
    if isempty(user_entry)
        user_entry = 'N';
    end
end

if strcmpi(user_entry,'Y') 
    try
        savepath
    catch ME
        warning('MFEToolbox:Path',['There was a problem saving the path.  The error was: \n' ME.message])
    end
end

if strcmpi(computer,'PCWIN64')
    user_entry = '';
    while ~ismember({user_entry},{'Y','N'}) && ~silent
        user_entry = input('Add MEX files to path? [Y(default)/N] ','s');
        user_entry = upper(user_entry);
        if isempty(user_entry)
            user_entry = 'Y';
        end
    end
    if strcmp(user_entry,'Y') || silent
        addpath(fullfile(pwd,mexDir));
    end
end

