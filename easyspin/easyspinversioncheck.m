% easyspinversioncheck Checks online for a new version 

% Possible calls:
% 
%   easyspinversioncheck
%   easyspinversioncheck(InstallInfo)
%   easyspinversioncheck(InstallInfo, Opt)
% 
%   [UpdateAvailable, OnlineVersion] = easyspinversioncheckchecks
%   [UpdateAvailable, OnlineVersion] = easyspinversioncheck(InstallInfo) 
%   [UpdateAvailable, OnlineVersion] = easyspinversioncheck(InstallInfo, Opt)
% 
% Input: 
%   InstallInfo             structure, get it by running: 
%                           InstallInfo = easyspininfo
%         .ReleaseChannel   release channel of current installation
%         .Version          EasySpin version
%         .Path             Installation path
% 
%   Opt                     options
%         .Silent           T/F - if T supress messages
%         .Branch           string, find newest version of specified
%                           release branch online
% 
% Output:
%   UpdateAvailable         T/F, T if update is available
%   OnlineVersion           string, most upto-date EasySpin version that is
%                           currently available for installed branch or
%                           branch that was specified with Opt.Branch

function varargout = easyspinversioncheck(VersionInfo, Opt)

messages = true;
testing = false;

if nargin == 2
  if isfield(Opt,'Silent') && Opt.Silent
    messages = false;
  end
  if isfield(Opt,'Test') && Opt.Test
    testing = true;
  end
end

% Check if MATLAB has write access to working directory
%--------------------------------------------------------------
fileName = '.\tmpeasyspinupdatetest.txt';
[fid,errmsg] = fopen(fileName, 'w');
if ~isempty(errmsg)
  if messages
    error('MATLAB can not write to the working directory. Please change your working directory and try again.');
  else
    fprintf('\n  MATLAB does not have write access to the current working directory.')
    fprintf('\n  For optimal performance, please switch to writeable working directory. \n');
    varargout = {0 []};
    return
  end
else
    fclose(fid);
    delete(fileName);
end

% First check if computer is online and easyspin.org reachable
%--------------------------------------------------------------
if ~testing
  if ispc
    [isOffline,~] = system('ping -n 1 www.google.com');
    [EasySpinOrgOffline,~] = system('ping -n 1 easyspin.org');
  elseif isunix
    [isOffline,~] = system('ping -c 1 www.google.com');
    [EasySpinOrgOffline,~] = system('ping -c 1 easyspin.org');
  end
  
  
  if isOffline
    msg = 'You need to be connect to the internet to check for a new EasySpin version or to update.';
    if messages; disp(msg); end
    varargout = {0 []};
    return
  elseif EasySpinOrgOffline
    msg = '<a href="easyspin.org">easyspin.org</a> appears to be offline, please try again later.';
    if messages; disp(msg); end
    varargout = {0 []};
    return
  end
end

% Computer is online, continue by processing info of installed version
%--------------------------------------------------------------
if nargin == 0
  % get version from easyspininfo if called without input argument
  VersionInfo = easyspininfo;
end

if nargin == 2 && isfield(Opt,'Branch')
  ReleaseChannel = Opt.Branch;
  InstalledVersion =  '0.0.0';
else
  ReleaseChannel   =  VersionInfo.ReleaseChannel;
  InstalledVersion =  VersionInfo.Version;
end

% stop looking for an update if EasySpin is on source control
if strcmp(ReleaseChannel,'$ReleaseChannel$')
  msg = 'Your EasySpin installation appears to be on version control, use your source control tool to check for a newer version.';
  if messages; disp(msg); end
  if nargout > 0
    varargout = {false []};
  end
  return
end

matchVersioning = '(\d+).(\d+).(\d+)-?([a-z]+)?.?(\d+)?';

VersionInfo = regexp(InstalledVersion,matchVersioning,'tokens');

if isempty(VersionInfo)
  error('Error polling currently installed version')
else
  VersionInfo = VersionInfo{:};
end

InstalledVersion = convertVersionNumber(VersionInfo);

% Get online versions and put them into a structure
%--------------------------------------------------------------
msg = 'Checking for a new version...';
if messages; disp(msg); end

if ~testing
  VersionsFile = websave('versions.txt','http://easyspin.org/easyspin/versions.txt');
else
  VersionsFile = './data/versions.txt';
end

fileID = fopen(VersionsFile,'r');

expression = ['([a-z]+):' matchVersioning];
tline = fgetl(fileID);

while ischar(tline)
  V = regexp(tline,expression,'tokens');
  if ~isempty(V)
    V = V{:};
    [OnlineVersions.(V{1}).Numeric, OnlineVersions.(V{1}).String] = convertVersionNumber(V(2:end));
  end
  tline = fgetl(fileID);
end

fclose(fileID);

if ~testing
  delete(VersionsFile);
end

% Check if the releasechannel of the currently installed version (or the
% requested releasechannel) are tracked on easyspin.org
%----------------------------------------------------------
if ~isfield(OnlineVersions,ReleaseChannel)
  msg = ['The release channel you are currently on - "' ReleaseChannel '" - does not support automatic updates.'];
  if messages; disp(msg); end
    varargout = {false []};
  return
end

% Compare versions if a new version is available
%----------------------------------------------------------
OnlineVersion = OnlineVersions.(ReleaseChannel).Numeric;

% Subtract local version from online version. Any version numbers of the
% local installation that are larger than the online versions, are
% negative. If online and local version are idential, Diffs is a vector of
% zeros.
Diffs = OnlineVersion - InstalledVersion;

% if both versions arent identical, Diffs contains nonzero elements
if any(Diffs > 0) && any(Diffs < 0)
  % if Diffs contains both negative and positive values (and zeros) their
  % positions must be compared. if the positive value has a lower index 
  % than the negative, a newer version is available to download
  UpdateAvailable = find(Diffs > 0,1) < find(Diffs < 0,1);
elseif any(Diffs > 0) && ~any(Diffs < 0)
  % if Diffs contains only positive values (and zeros), a newer version is
  % available
  UpdateAvailable = true;
elseif ~any(Diffs > 0) && any(Diffs < 0)
  % if Diffs contains only negative values (and zeros), the local version
  % is up to date
  UpdateAvailable = false;
else
  % if both versions are identical, all elements in Diffs are zero
  UpdateAvailable = false;
end

if UpdateAvailable
  msg = ['A new EasySpin version (' OnlineVersions.(ReleaseChannel).String ') is available online.' newline];
  msg = [msg 'Type and run "easyspinupdate" to update.'];
else
  msg = 'You are using the most recent version of EasySpin.';
end
if messages; disp(msg); end

if nargout > 0
  varargout = {UpdateAvailable OnlineVersions.(ReleaseChannel).String};
end

function [VersionArray, VersionString] = convertVersionNumber(VersionInfo)
% Translates the cell array, obtained from regexp, to a numeric value
% that encodes the EasySpin version
% This version number, as well as the semantic version number in string
% form is returned

VersionString = join(VersionInfo(1:3),'.');

if isempty(VersionInfo{4})
  VersionInfo{4} = '5';
else
  VersionString = join({VersionString{:} VersionInfo{4}},'-');
  switch VersionInfo{4}
    case 'dev'
      VersionInfo{4} = '0';
    case 'alpha'
      VersionInfo{4} = '1';
    case 'beta'
      VersionInfo{4} = '2';
  end
end

if isempty(VersionInfo{5})
  VersionInfo{5} = '0';
else
  VersionString = join({VersionString{:} VersionInfo{5}},'.');
end

VersionArray = cellfun(@str2double,VersionInfo);

VersionString = VersionString{1};

