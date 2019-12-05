function esbuild
%========================================================================
%                 EasySpin build script
%========================================================================


%========================================================================
% Build settings
%========================================================================
% ReleaseID: MAJOR.MINOR.PATCH
%   MAJOR: Change only when major new functionality is implemented
%     including incompatible changes.
%   MINOR: Change when new functionality is added in a backwards-
%     compatible manner.
%   PATCH: Increment for every bugfix release.
% Roughly follow guidelines of semantic versioning, see http://semver.org/
ReleaseID = '%ReleaseID%'; % major.minor.patch
ReleaseChannel = '%ReleaseChannel%'; % release channel, used for automatic updates

% Expiry date of release, see eschecker.m
% Months to add:
MonthsToExpiry = %Months; % Will be replaced by build.pl

[Year, Month, Day] = datevec(now);
ExpiryDate = datestr(datenum(Year,Month+MonthsToExpiry+1,1)-1);

% Cutoff date for date checking, see eschecker.m
HorizonDate = datestr(datenum(Year+4,12,31));

% Folders
SourceDir = ['%SourceDir%'];
ZipDestDir = ['%ZipDestDir%'];

%========================================================================


%------------------------------------------------------------------------
clc
v = sscanf(version,'%f',1);
% if v>8.4
%   error('EasySpin build must be done with Matlab 8.4 (R2014b).');
% end

%error('Must include perl script that replaces $ReleaseID$ and $ReleaseDate$ globally.');

BuildTimeStamp = datestr(now,'yyyymmdd-HHMMSS');
% BuildID = sprintf('%s+%s',ReleaseID,BuildTimeStamp);
BuildID = sprintf('%s',ReleaseID);
ReleaseDate = sprintf('%04d-%02d-%02d',Year,Month,Day);

fprintf('Building EasySpin %s.\n',BuildID);

fprintf('Checking for source folders...');

Dirs = {'docs','easyspin','examples'};
for k = 1:numel(Dirs)
  if ~exist([SourceDir filesep Dirs{k}],'dir')
    error('Could not find EasySpin source subdirectories.');
  end
end

fprintf(' ok\n');

fprintf('Creating build folder...');

BuildFolder = [tempdir 'easyspin-' ReleaseID];

if exist(BuildFolder,'dir')
  ok = rmdir(BuildFolder);
  if (~ok)
    error('Could not remove old build folder %s.',BuildFolder);
  end
end
if exist(BuildFolder,'file')
  delete(BuildFolder);
end

ok = mkdir(BuildFolder);
if (~ok)
  error('Could not create build folder %s.',BuildFolder);
end

fprintf(' ok\n');
fprintf('  Temp folder: %s\n',tempdir);
fprintf('  Build folder: %s\n',BuildFolder);

TbxFolder = [BuildFolder filesep 'easyspin'];
ExmplFolder = [BuildFolder filesep 'examples'];
DocFolder = [BuildFolder filesep 'documentation'];

TbxPcodeDir = [BuildFolder filesep 'easyspinpcode'];

%------------------------------------------------------------------------
% Toolbox folder
%------------------------------------------------------------------------
disp('Toolbox');
TbxSrcDir = [SourceDir filesep 'easyspin'];

% Assert that .c files do not contain // comments (which are not
% supported by strict ANSI-C compilers, but are supported by the C99
% standard .)
enforceANSICcomments = false;
if enforceANSICcomments
  fprintf('  Checking *.c files for absence of // comments...');
  cfiles = dir([TbxSrcDir filesep 'private' filesep '*.c']);
  for k=1:numel(cfiles)
    Lines = textread([TbxSrcDir filesep 'private' filesep cfiles(k).name],'%s','whitespace','\n');
    f = [];
    for q=1:numel(Lines)
      f = [f strfind(Lines{q},'//')];
    end
    if ~isempty(f)
      error('Found // comment in file %s',cfiles(k).name);
    end
  end
  fprintf(' ok\n');
end

fprintf('  Generating toolbox folder...');
mkdir(BuildFolder,'easyspin');
fprintf(' ok\n');

fprintf('  Copying toolbox files... ');
copyfile(TbxSrcDir,TbxFolder);
fprintf(' ok\n');

% Release tag and expiry date
fprintf('  Setting release information and expiry date... ');
replacestr([TbxFolder filesep 'info.xml'],'$ReleaseID$',ReleaseID);
replacestr([TbxFolder filesep 'easyspininfo.m'],'$ReleaseID$',ReleaseID);
replacestr([TbxFolder filesep 'easyspininfo.m'],'$ReleaseChannel$',ReleaseChannel);
replacestr([TbxFolder filesep 'easyspininfo.m'],'$ReleaseDate$',ReleaseDate);
replacestr([TbxFolder filesep 'easyspininfo.m'],'$ExpiryDate$',ExpiryDate);
replacestr([TbxFolder filesep 'eschecker.m'],'888888',num2str(datenum(ExpiryDate)));
replacestr([TbxFolder filesep 'eschecker.m'],'999999',num2str(datenum(HorizonDate)));
fprintf(' ok\n');


%------------------------------------------------------------------------
% P-coding
%------------------------------------------------------------------------
% Copy everything to a third directory for pcoding
% This makes sure .m files are not newer than .p files
%  (otherwise Matlab complains about potentially obsolete p-files)
fprintf('  Copying to a new folder for p-coding...');
mkdir(BuildFolder,'easyspinpcode');
copyfile(TbxFolder,TbxPcodeDir);
fprintf(' ok\n');

fprintf('  Extracting help from all m-files...');
extracthelp(TbxFolder);
fprintf(' ok\n');

fprintf('  Deleting m-files in private subdirectory...');
delete([TbxFolder filesep 'private' filesep '*.m']);
fprintf(' ok\n');

fprintf('  P-coding all m-files...');
curdir = pwd;
cd(TbxFolder);
pcode([TbxPcodeDir filesep '*.m']);
cd([TbxFolder filesep 'private']);
pcode([TbxPcodeDir filesep 'private' filesep '*.m']);
cd(curdir);
fprintf(' ok\n');

rmdir(TbxPcodeDir,'s');


%------------------------------------------------------------------------
% Examples
%------------------------------------------------------------------------
disp('Examples');

fprintf('  generating examples dir and copying files...');
mkdir(BuildFolder,'examples');
copyfile([SourceDir filesep 'examples'],ExmplFolder);
fprintf(' ok\n');

%------------------------------------------------------------------------
% Documentation
%------------------------------------------------------------------------
disp('Documentation');

fprintf('  generating documentation folder...');
mkdir(BuildFolder,'documentation');
fprintf(' ok\n');

fprintf('  copying files...');
DocSrcDir = [SourceDir filesep 'docs'];
copyfile(DocSrcDir,DocFolder);
fprintf(' ok\n');

% Replace all $ReleaseID$ in html files with actual release ID
fprintf('  inserting release ID into documentation...');
docFiles = dir([DocFolder filesep '*.html']);
for iFile = 1:numel(docFiles)
  replacestr([DocFolder filesep docFiles(iFile).name],'$ReleaseID$',ReleaseID);
end
fprintf(' ok\n');


%------------------------------------------------------------------------
% Packaging
%------------------------------------------------------------------------
disp('Packaging');
% package for public release
if ZipDestDir(end)==filesep, ZipDestDir(end) = []; end
% zipFile = [ZipDestDir filesep 'easyspin-' BuildID '.zip'];
zipFileShort = [ZipDestDir filesep 'easyspin-' ReleaseID '.zip'];

% fprintf('  packing zip file %s...',zipFile);
% if exist(zipFile,'file')
%   try
%     delete(zipFile);
%   catch
%     error('Cannot delete zip file.');
%   end
% end
% zip(zipFile,BuildFolder);
% fprintf(' ok\n');

fprintf('  packing zip file %s...',zipFileShort);
zip(zipFileShort,BuildFolder);
fprintf(' ok\n');


%------------------------------------------------------------------------
% Clean-up
%------------------------------------------------------------------------
fprintf('Removing build tree...');
rmdir(BuildFolder,'s');
fprintf(' ok\n');

disp('Done!');

return


%========================================================================
%========================================================================
%========================================================================
%========================================================================

% extracthelp  Extracts help from m files
%
%    extracthelp(myDir)
%
%    Extract the help from all m-files in the directory
%    myDir, i.e. removes all the code.

function extracthelp(pDir)

mFiles = dir([pDir filesep '*.m']);
mFiles = sort({mFiles.name});
for k = 1:length(mFiles)
  mFiles{k} = [pDir filesep mFiles{k}];
end

% Copy help lines from dir/file.m to dir2/file.m
for k = 1:length(mFiles)
  %disp(['  help extraction of ' mFiles{k}]);

  infile = fopen(mFiles{k},'r');
  CommentTag = '%';
  ok = 1;
  LineStr = {};
  while ok
    thisLine = fgetl(infile);
    if isempty(thisLine) || all(thisLine==-1); break; end
    ok = thisLine(1)==CommentTag;
    if ok, LineStr{end+1} = thisLine; end
  end
  fclose(infile);

  outfile = fopen(mFiles{k},'w');
  for q=1:numel(LineStr)
    fprintf(outfile,'%s\n',LineStr{q});
  end
  fclose(outfile);
end


function replacestr(fname,S0,S1)

fid = fopen(fname);
Lines = textscan(fid,'%s','whitespace','','delimiter','\n');
Lines = Lines{1};
fclose(fid);
nLines = numel(Lines);

n = length(S0);
found = false;
for k = 1:nLines
  q = strfind(Lines{k},S0);
  if ~isempty(q)
    Lines{k} = [Lines{k}(1:q-1) S1 Lines{k}(q+n:end)];
    found = true;
  end
end

if found
  fid = fopen(fname,'w');
  for k=1:nLines
    fprintf(fid,'%s\n',Lines{k});
  end
  fclose(fid);
end
