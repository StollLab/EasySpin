% easyspin_info  Installation information and test 
%
%   Tests the EasySpin installation and prints details.

function varargout = easyspin_info()

Display = nargout==0;

warning(chkmlver);

MATLABversion = builtin('version');

% Determine EasySpin paths
%-------------------------------------------------------------------------------
esPath = fileparts(which(mfilename));
esRoot = esPath(1:end-length('\easyspin\private'));
esFunctionsFolder = [esRoot filesep 'easyspin'];
allFiles = what(esFunctionsFolder);


% Get release information
%-------------------------------------------------------------------------------
esVersion = '$ReleaseID$';  % replaced with actual version by build script
esDate = '$ReleaseDate$'; % replaced with actual date by build script
esExpiryDate = '$ExpiryDate$'; % replaced with actual date by build script
esReleaseChannel = '$ReleaseChannel$'; % replaced with release channel by build script


% Determine operating system
%-------------------------------------------------------------------------------
if isunix
  [~,platform] = unix('uname -s -r');
else
  [~,platform] = dos('ver');
end
platform(platform<31) = [];


% Prepare output
%-------------------------------------------------------------------------------
VersionInfo.ReleaseChannel = esReleaseChannel;
VersionInfo.Version = esVersion;
VersionInfo.Root = esRoot;
VersionInfo.Date = esDate;
VersionInfo.ExpiryDate = esExpiryDate;
VersionInfo.MATLABversion = MATLABversion;
VersionInfo.Platform = platform;
if nargout>0
  varargout = {VersionInfo};
  return
end


% Check if mex files are available, compile if necessary
%-------------------------------------------------------------------------------
% Get list of mex files from C source files in private subfolder
SrcFiles = dir([esFunctionsFolder filesep 'private' filesep '*.c']);
nSrcFiles = numel(SrcFiles);
for k = 1:nSrcFiles
  mexFiles{k} = SrcFiles(k).name(1:end-2);
end

% Check whether mex files are present
for k = 1:nSrcFiles
  mexed(k) = exist(mexFiles{k})==3;
end

% If not, try to compile and check again
if ~all(mexed)
  fprintf('\nEasySpin not yet compiled for this platform (%d/%d files). Attempting to compile.\n\n',sum(mexed),numel(mexed));
  easyspin_compile;
  for k = 1:nSrcFiles
    mexed(k) = exist(mexFiles{k})==3;
  end
  if all(mexed)
    fprintf('\nEasySpin is now compiled for this platform.\n\n');
  else
    fprintf('\nEasySpin compilation unsuccessful.\n\n');
  end
end
clear functions


% Display information
%-------------------------------------------------------------------------------
if Display
  fprintf('==================================================================\n');
  fprintf('Information about the installed EasySpin version\n');
  fprintf('==================================================================\n');
  isDevelopmentVersion = esVersion(1)=='$';
  if isDevelopmentVersion
    fprintf('  Version:          development\n');
    fprintf('  Expiration date:  none\n');
  else
    fprintf('  Version:          %s (%s)\n',esVersion,esDate);
    fprintf('  Expiration date:  %s\n',esExpiryDate);
  end
end

if Display
  fprintf('  Root folder:      %s\n',esRoot);
  fprintf('  mex files:        ');
  if all(mexed)
    fprintf([mexext, ', all present']);
  else
    fprintf([mexext, ', %d/%d missing'],sum(~mexed)/numel(mexed));
  end
  fprintf('\n');

  % Display information about MATLAB and system
  fprintf('------------------------------------------------------------------\n');
  fprintf('  MATLAB version:   %s\n',VersionInfo.MATLABversion);
  fprintf('  Platform:         %s\n',platform);
  fprintf('  System date:      %s\n',char(datetime));
  fprintf('  Temp dir:         %s\n',tempdir);
  fprintf('------------------------------------------------------------------\n');

  % Display information about how to access the documentation
  docEntry = easyspin_doc;
  fprintf('  To access the documentation:\n');
  fprintf('  - Click <a href="%s">here</a>\n',docEntry);
  fprintf('  - In the command window, type easyspin doc\n')
  fprintf('  - Point your browser to <a href="%s">%s</a>\n',docEntry,docEntry);
  fprintf('==================================================================\n');
  
  % Check online for update
  %-----------------------------------------------------------------------------
  UpdateOpt.Silent = true;
  [UpdateAvailable, NewerVersion] = easyspinversioncheck(VersionInfo,UpdateOpt);
  if UpdateAvailable
    msg = [newline '   A new EasySpin version (' NewerVersion ') is available online.' newline];
    msg = [msg 'Type and run "easyspinupdate" to update.'];
    disp(msg)
  end
  
end


% Determine whether EasySpin's directory is on the search path
%-------------------------------------------------------------------------------
onSearchPath = contains(upper([path pathsep]),upper([esFunctionsFolder pathsep]));
if ~onSearchPath
  if Display
    fprintf('\n  The EasySpin folder %s is not in MATLAB''s search path. Please add it.\n',esFunctionFolder);
  end
end


% Check for name conflicts
%-------------------------------------------------------------------------------
if Display
  %fprintf('    Checking for name conflicts... ');
end
Shadowing = {};
Shadowed = {};
for iFunction = 1:length(allFiles.m)
  if ~isempty(strfind(allFiles.m{iFunction},'Contents')); continue; end
  Instances = which(allFiles.m{iFunction},'-all');

  % Remove MATLAB class methods
  % Class methods are only called for objects of their class.
  % Private functions are only called from functions in the parent folder.
  classMethodID = [filesep '@'];
  privateFunctionID = [filesep 'private' filesep];
  for k = numel(Instances):-1:1
    isClassMethod = contains(Instances{k},classMethodID);
    isPrivateFunction = contains(Instances{k},privateFunctionID);
    if isClassMethod || isPrivateFunction
      Instances(k) = [];
    end
  end
  
  if numel(Instances)>1
    FirstPath = fileparts(Instances{1});
    if strcmp(FirstPath,esFunctionsFolder)
      Shadowed = [Shadowed; Instances(2:end)];
    else
      Shadowing = [Shadowing; Instances(1)];
    end
  end
  
end


% List files that are shadowed by EasySpin
%-------------------------------------------------------------------------------
listShadowed = false; % don't list ...
if listShadowed && ~isempty(Shadowed)
  fprintf('\n  EasySpin functions shadow the following functions\n');
  for k=1:length(Shadowed)
    fprintf('    > %s\n',Shadowed{k});
  end
  fprintf('  If these functions should work, either move or rename them.');
  fprintf('\n  Otherwise remove them from MATLAB''s search path or remove them.');
end


% List files that shadow EasySpin's functionality
%-------------------------------------------------------------------------------
if ~isempty(Shadowing)
  fprintf('\n  The following functions shadow EasySpin functions with the same name\n');
  for k=1:length(Shadowing)
    fprintf('    > %s\n',Shadowing{k});
  end
  fprintf('\n  EasySpin works properly only if you rename, move or remove these\n');
  fprintf('  functions. Alternatively, remove the corresponding directory from\n');
  fprintf('  the MATLAB path.\n');
end

if isempty(Shadowed) && isempty(Shadowing)
  if Display
    %fprintf('ok\n');
  end
end

end
