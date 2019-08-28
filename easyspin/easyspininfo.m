% easyspininfo  Installation information and test 
%
%   Tests the EasySpin installation and prints
%   release details.

function varargout = easyspininfo()

Display = (nargout==0);

error(chkmlver);

% Determine EasySpin path.
%--------------------------------------------------------------
esPath = fileparts(which(mfilename));
AllFiles = what(esPath);

% Get release information.
%--------------------------------------------------------------
esVersion = '$ReleaseID$';  % is replaced with actual version by build script
esDate = '$ReleaseDate$'; % is replaced with actual date by build script
esExpiryDate = '$ExpiryDate$'; % is replaced with actual date by build script
esReleaseChannel = '$ReleaseChannel$'; % is replaced with release channel by build script



out.ReleaseChannel = esReleaseChannel;
out.Version = esVersion;
out.Path = esPath;
if (nargout>0)
  varargout = {out};
end


% Determine operating system
%--------------------------------------------------------------
if isunix
  [ignore,platform] = unix('uname -s -r');
else
  [ignore,platform] = dos('ver');
end
platform(platform<31) = [];

% Check if mex files are available, compile if necessary
%-------------------------------------------------------------
% Get list of mex files from C source files in private subfolder
SrcFiles = dir([esPath '/private/*.c']);
for k = 1:numel(SrcFiles)
  MexFiles{k} = SrcFiles(k).name(1:end-2);
end

% Check whether mex files are present
for k = 1:numel(SrcFiles)  
  MexFilesPresent(k) = exist(MexFiles{k})==3;
end

% If not, try to compile and check again
if ~all(MexFilesPresent)
  fprintf('\nEasySpin not yet compiled for this platform (%d/%d files). Attempting to compile.\n\n',sum(MexFilesPresent),numel(MexFilesPresent));
  easyspincompile;
  for k = 1:numel(MexFiles)
    MexFilesPresent(k) = exist(MexFiles{k})==3;
  end
  if all(MexFilesPresent)
    fprintf('\nEasySpin is now compiled for this platform.\n\n');
  else
    fprintf('\nEasySpin compilation unsuccessful.\n\n');
  end
end
clear functions

% Display information
%--------------------------------------------------------------
if Display
  fprintf('==================================================================\n');
  fprintf('  Release:          %s (%s)\n',esVersion,esDate);
end

Diagnostics = true;
if Diagnostics && Display
  fprintf('  Expiry date:      %s\n',esExpiryDate);
  fprintf('  Folder:           %s\n',esPath);
  fprintf('  MATLAB version:   %s\n',builtin('version'));
  fprintf('  Platform:         %s\n',platform);
  MexFiles = dir([esPath filesep 'private' filesep '*.c']);
  for k=1:numel(MexFiles)
    mexed(k) = exist(MexFiles(k).name(1:end-2))==3;
  end
  fprintf('  mex-files:        ');
  if all(mexed)
    fprintf([mexext, ', ok']);
  else
    fprintf([mexext, ', missing']);
  end
  fprintf('\n');
  fprintf('  System date:      %s\n',datestr(now));
  fprintf('  Temp dir:         %s\n',tempdir);
  fprintf('==================================================================\n');
  
  % Check online for update
  %----------------------------------------------------------
  if ispc
    [isOffline,~] = system('ping -n 1 www.google.com');
    [EasySpinOrgOffline,~] = system('ping -n 1 easyspin.org');
  elseif isunix
    [isOffline,~] = system('ping -c 1 www.google.com');
    [EasySpinOrgOffline,~] = system('ping -c 1 easyspin.org');
  end
  
  isOnline = ~isOffline;
  EasySpinOrgOnline = ~EasySpinOrgOffline;
  
  if isOnline && EasySpinOrgOnline
    easyspinupdate(out);
  end
end


% Determine whether EasySpin's directory is on the search path.
%--------------------------------------------------------------
OnSearchPath = strfind(upper([path pathsep]),upper([esPath pathsep]));

if isempty(OnSearchPath)
  if Display
    fprintf('The EasySpin folder is not in MATLAB''s search path. Please add it.\n\n');
  end
end

% Check for name conflicts.
%--------------------------------------------------------------
if Display
  %fprintf('    Checking for name conflicts... ');
end
Shadowing = {};
Shadowed = {};
for iFunction = 1:length(AllFiles.m)
  if ~isempty(strfind(AllFiles.m{iFunction},'Contents')); continue; end
  Instances = which(AllFiles.m{iFunction},'-all');

  % Remove MATLAB class methods. They are only called for objects
  % of their class.
  for k = numel(Instances):-1:1
    if strfind(Instances{k},[filesep '@'])
        Instances(k) = [];
    end
  end

  if numel(Instances)>1
    FirstPath = fileparts(Instances{1});
    if strcmp(FirstPath,esPath)
      Shadowed = [Shadowed; Instances(2:end)];
    else
      Shadowing = [Shadowing; Instances(1)];
    end
  end
end

% List files that are shadowed by EasySpin.
%--------------------------------------------------------------
Shadowed = {}; % don't list ...
if ~isempty(Shadowed)
  fprintf('\n\n    EasySpin functions shadow the following functions\n');
  for k=1:length(Shadowed)
    fprintf('    > %s\n',Shadowed{k});
  end
  fprintf('    If these functions should work, either move or rename them.');
  fprintf('\n    Otherwise remove them from MATLAB''s search path or remove them.');
end

% List files that shadow EasySpin's functionality.
%--------------------------------------------------------------
if ~isempty(Shadowing)
  fprintf('\n   The following functions shadow EasySpin functions with the same name\n');
  for k=1:length(Shadowing)
    fprintf('    > %s\n',Shadowing{k});
  end
  fprintf('\n    EasySpin works properly only if you rename, move or remove these\n');
  fprintf('    functions. Alternatively, remove the corresponding directory from\n');
  fprintf('    the MATLAB path.\n\n');
end

if isempty(Shadowed) & isempty(Shadowing)
  if Display
    %fprintf('ok\n');
  end
else
end

return
