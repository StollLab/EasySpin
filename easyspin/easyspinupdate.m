% esupdate  Checks online for a new version and updateds if requested 
%
% Possible calls:
%   esupdate           checks for update and installs update
%   esupdate('check')  only checks if there is a newer version available
%   esupdate('stable') updates to the most recent version for the specified
%                      release channel, available options are: 'stable',
%                      'legacy', 'development'
%   esupdate('5.2.22') installs the specified version of EasySpin, if the
%                      version is not found online, an error message is
%                      returned

function varargout = easyspinupdate(varargin)

% release channels that are kept tracked of on easyspin.org
AllowedReleaseChannels = {'stable','development','legacy'};

ForceUpdate = false;
CheckVersion = false;
UpdateToVersion = false;
NewVersionAvailable = false;

% Get information from current install
%--------------------------------------------------------------
if nargin == 1 && isstruct(varargin{1})
  % easyspininfo calls esupdate with required information already
  installedES = varargin{1};
else
  % obtains information from easyspininfo;
  installedES = easyspininfo;
end

ReleaseChannel = installedES.ReleaseChannel;
InstalledVersion =  installedES.Version;
InstallationPath = installedES.Path;

% Get the real installation path, which is two directories above the easyspin functions:
Path = strsplit(InstallationPath,filesep);
InstallationPath = join(Path(1:end-2),filesep);
InstallationPath = InstallationPath{1};

% In case easyspininfo was not correctly tagged during the build process
if strcmp(ReleaseChannel,'$ReleaseChannel$')
  ReleaseChannel = 'stable';
end

if strcmp(InstalledVersion,'$ReleaseID$')
  InstalledVersion = '0.0.0';
end


% Process input arguments if any
%---------------------------------------------------------- 
if nargin > 0 
  
  if isstruct(varargin{1}) || strcmp(varargin{1},'check')
    % only check for an update
    CheckVersion = true;
  else
    Call = varargin{1};    
    ForceUpdate = true;
    
    if any(strcmp(Call,AllowedReleaseChannels))
      % input corresponds to a release channel
      ReleaseChannel = Call;
    else
      % input is a version
      UpdateToVersion = Call;
    end
  end
  
end

if ~UpdateToVersion
  % Translate currently installed version to a numeric value
  %----------------------------------------------------------
  matchVersioning = '(\d+).(\d+).(\d+)-?([a-z]+)?.?(\d+)?';

  VersionInfo = regexp(InstalledVersion,matchVersioning,'tokens');

  if isempty(VersionInfo)
    error('Error polling currently installed version')
  else
    VersionInfo = VersionInfo{:};
  end

  installedVersionNumber = getVersionNumber(VersionInfo);

  % Get versions that are available online
  %----------------------------------------------------------
  msg = 'Checking for a new version online';
  disp(msg)
  VersionsFile = websave('versions.txt','http://easyspin.org/easyspin/versions.txt');
  fileID = fopen(VersionsFile,'r');
  expression = ['([a-z]+):' matchVersioning];

  tline = fgetl(fileID);
  while ischar(tline)
      V = regexp(tline,expression,'tokens');
      if ~isempty(V)
        V = V{:};
        [OnlineVersions.(V{1}).Numeric, OnlineVersions.(V{1}).String] = getVersionNumber(V(2:end));
      end
      tline = fgetl(fileID);
  end

  fclose(fileID);
  delete(VersionsFile);

  % Check if the releasechannel of the currently installed version (or the
  % requested releasechannel) are tracked on easyspin.org
  %----------------------------------------------------------
  if ~isfield(OnlineVersions,ReleaseChannel)
    errMsg = ['The release channel "' ReleaseChannel '" was not found on easyspin.org.'];
    error(errMsg);
  end

  % See if a new version is available
  %----------------------------------------------------------
  NewVersionAvailable = OnlineVersions.(ReleaseChannel).Numeric > installedVersionNumber;
  
  if CheckVersion
    varargout = {NewVersionAvailable};
    msg = ['A new EasySpin version (' OnlineVersions.(ReleaseChannel).String ') is available online.' newline];
    msg = [msg 'Type and run "esupdate" to update.'];
    disp(msg);
    % Terminate if esupdate was only called to check the version
    return
  end
  
end

% Download and install new version
%----------------------------------------------------------
if NewVersionAvailable || ForceUpdate
  
  % Decide which version to download
  if UpdateToVersion
    VersionToGet = UpdateToVersion;
  else
    VersionToGet = OnlineVersions.(ReleaseChannel).String;
  end
  
  disp(['Downloading new EasySpin version (' VersionToGet ')']);
  zipName = ['easyspin-' VersionToGet '.zip'];
  
  % download from easyspin.org
  try
    zipFile = websave(zipName,['http://easyspin.org/easyspin/test/' zipName]);
  catch
    delete([zipName '.html']); % if the file can not be downloaded, MATLAB 
                               % creates a file 'filename.html', this
                               % removes the file
    errMsg = ['The file ' zipName ' was not found on easyspin.org.'];
    error(errMsg);
  end
  
  disp('Installing...');
  
  % unzip to destination
  Destination = [InstallationPath filesep];
  unzip(zipFile,Destination);
  
  % remove downloaded zip
  delete(zipFile);
  
  % Add to Path
  NewESPath = [Destination 'easyspin-' VersionToGet filesep 'easyspin' filesep];
  
  % If for whatever reason the path contains an error and is not available
  % on the system, an error is returned, notifiying the user that EasySpin
  % was downloaded, but could not be added to the path
  if isfolder(NewESPath)
    addpath(NewESPath);
    savepath
    msg = ['EasySpin was succesfully installed to ' newline NewESPath newline 'and added to Paths.'];
    disp(msg);
  else
    errMsg = ['EasySpin was succecsfully downloaded to ' newline NewESPath newline];
    errMsg = [errMsg 'But adding it to the path failed. Please do so manually.'];
    error(errMsg)
  end
else
  disp('You are using the most recent version of EasySpin.')
end


function [VersionNumber, VersionString] = getVersionNumber(VersionInfo)
  % Translates the cell array, obtained from regexp, to a numeric value
  % that encodes the EasySpin version
  % This version number, as well as the semantic version number in string
  % form is returned
  
  VersionNumber = 10000*str2double(VersionInfo{1})+100*str2double(VersionInfo{2})+str2double(VersionInfo{3});
  VersionString = join(VersionInfo(1:3),'.');
  
  % development version:
  if length(VersionInfo) > 3 && ~isempty(VersionInfo{4})
    % this takes care of the alpha, beta and dev versions
    switch VersionInfo{4}
      case 'alpha'
        VersionNumber = VersionNumber + 0.2;
      case 'beta'
        VersionNumber = VersionNumber + 0.3;
      case 'dev'
        VersionNumber = VersionNumber + 0.1;
    end
    VersionString = join({VersionString{:} VersionInfo{4}},'-');
    
    % in case the development version ends with a numeric value, e.g.
    % 6.0.0-dev.3 vs 6.0.0-dev.4
    if length(VersionInfo) > 3 && ~isempty(VersionInfo{5})
      VersionNumber = VersionNumber + 0.01*str2double(VersionInfo{5});
      VersionString = join({VersionString{:} VersionInfo{5}},'.');
    end
    
  end
  VersionString = VersionString{1};
  