% easyspinupdate  Checks online for a new version and updateds if requested 
%
% Possible calls:
%   easyspinupdate           checks online for update and downloads and
%                            install newest version
%   easyspinupdate('stable') downloads most recent version of requested
%                            release branch (stable or development)
%   easyspinupdate('5.2.22') downloads specific EasySpin version

function varargout = easyspinupdate(OnlineVersion)

% ---------------------------------------------------------------
% If easyspinupdate is called without argument, a version check is run.
% If a newer version is available, easyspinupdate is called again with the
% newer version as argument.
if nargin == 0
  [UpdateAvailable, OnlineVersion] = easyspinversioncheck;
  if UpdateAvailable
    easyspinupdate(OnlineVersion);
  end
  return
end

if all(isstrprop(OnlineVersion,'alpha'))
  InstalledVersion = easyspininfo;
  UpdateOpt.Branch = OnlineVersion;
  UpdateOpt.Silent = true;
  [~, OnlineVersion] = easyspinversioncheck(InstalledVersion,UpdateOpt);
  if isempty(OnlineVersion)
    msg = [UpdateOpt.Branch ' is not a valid branch name.'];
    disp(msg)
    return
  end
end

% ----------------------------------------------------------------
% First check if server can be reached
if ispc
  [isOffline,~] = system('ping -n 1 www.google.com');
  [EasySpinOrgOffline,~] = system('ping -n 1 easyspin.org');
elseif isunix
  [isOffline,~] = system('ping -c 1 www.google.com');
  [EasySpinOrgOffline,~] = system('ping -c 1 easyspin.org');
end

if isOffline
  msg = 'You have to be connect to the internet to update EasySpin.';
  disp(msg)
  return
end

if EasySpinOrgOffline
  msg = '<a href="easyspin.org">easyspin.org</a> appears to be offline, please try again later.';
  disp(msg)
  return
end

% ---------------------------------------------------------------
% Download and install

VersionToGet = OnlineVersion;

% Determine installation path of currently installed EasySpin
InstalledVersion = easyspininfo;
InstallationPath = InstalledVersion.Path;

% The installation target is two directories above the easyspin functions:
Path = strsplit(InstallationPath,filesep);
InstallationPath = join(Path(1:end-2),filesep);
InstallationPath = InstallationPath{1};

% Before downloading, check if the previous installation is in a write
% protected directory
fileName = join([InstallationPath "easyspininstalltest.txt"],filesep);
[fid,errmsg] = fopen(fileName, 'w');
if ~isempty(errmsg) 
    error(['The EasySpin installation directory (' InstallationPath ') appears to be write protected. Please move EasySpin to a different directory and retry update or manually install the new EasySpin version.']);
else
    fclose(fid);
    delete(fileName);
end

OldPath = join(Path(1:end-1),filesep);
OldPath = OldPath{1};

disp(['Downloading EasySpin version (' VersionToGet ')']);
zipName = ['easyspin-' VersionToGet '.zip'];

% set time out to 60 seconds
options = weboptions('Timeout',60);

% download from easyspin.org
try
  zipFile = websave(zipName,['http://easyspin.org/easyspin/' zipName],options);
catch
  if exist([zipName '.html'], 'file') == 2
    % if the file can not be downloaded, MATLAB
    % creates a file 'filename.html', this removes the file
    delete([zipName '.html']); 
    errMsg = ['The file ' zipName ' was not found on easyspin.org.'];
  else
    errMsg = ['It appears the connection has timed out. Please try again.'];
  end
  error(errMsg);
end

disp('Installing...');
  
% unzip to destination
Destination = [InstallationPath filesep];
unzip(zipFile,Destination);

% remove downloaded zip
delete(zipFile);

% ---------------------------------------------------------------
% Add to Path and clean up
newESPath = [Destination 'easyspin-' VersionToGet filesep 'easyspin' filesep];

if exist(newESPath,'dir')
  addpath(newESPath);
  savepath
  msg = ['EasySpin was succesfully installed to ' newline newESPath newline 'and added to the MATLAB search paths.' newline];
  msg = [msg 'For optimal perfomance, your should remove EasySpin installation (' OldPath ') from the MATLAB search paths and delete the folder from your system.']; 
  disp(msg);
else
  errMsg = ['EasySpin was succecsfully downloaded to ' newline newESPath newline];
  errMsg = [errMsg 'But adding it to the path failed. Please do so manually.'];
  error(errMsg);
end
