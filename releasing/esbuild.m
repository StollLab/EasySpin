function esbuild
%===========================================================
%                 EasySpin build script
%===========================================================

ReleaseID = '5.1.0'; % major.minor.patch
betaVersion = true;
packDevtree = false;
ExpiryDate = '30-Jun-2017';
HorizonDate = '31-Dec-2021';

baseDir = 'C:\Users\abc\Documents\work';
SourceDir = [baseDir '\easyspin-dev'];
ZipDestDir = baseDir;

%-----------------------------------------------------------
clc
v = sscanf(version,'%f',1);
if (v>7.5)
  error('EasySpin build must be done with Matlab 7.5 (R2007b, because of pcode).');
end

%error('Must include perl script that replaces $ReleaseID$ and $ReleaseDate$ globally.');

[Y,M,D,H,MI,S] = datevec(now);
ReleaseDate = sprintf('%04d-%02d-%02d',Y,M,D);
BuildTime = sprintf('%04d%02d%02d%02d%02d',Y,M,D,H,MI);
% Matlab 6.5 equivalent for
%BuildID = datestr(now,'yymmddHHMM');

% This file contains the last BuildNumber
load lastbuildid
BuildNumber = BuildNumber + 1;
if betaVersion
  BuildID = sprintf('%s-beta+%d',ReleaseID,BuildNumber);
else
  BuildID = sprintf('%s+%d',ReleaseID,BuildNumber);
end

fprintf('Building EasySpin %s.\n',BuildID);

disp('Checking for source directories.');

Dirs = {'docs','easyspin','examples'};
for k = 1:numel(Dirs)
  if ~exist([SourceDir filesep Dirs{k}],'dir')
    error('Could not find EasySpin source subdirectories.');
  end
end

disp('Creating build tree.');

if betaVersion
  BuildDir = [tempdir 'easyspin-' BuildID];
else
  BuildDir = [tempdir 'easyspin-' ReleaseID];
end
fprintf('Builddir is %s\n',BuildDir);
if exist(BuildDir,'dir');
  try
    rmdir(BuildDir,'s');
  catch
    error('Could not remove old build folder %s.',BuildDir);
  end
end
if betaVersion
  mkdir(tempdir,['easyspin-' BuildID]);
else
  mkdir(tempdir,['easyspin-' ReleaseID]);
end

TbxDir = [BuildDir filesep 'easyspin'];
ExmplDir = [BuildDir filesep 'examples'];
DocDir = [BuildDir filesep 'documentation'];

TbxPcodeDir = [BuildDir filesep 'easyspinpcode'];

%------------------------------------------------------------
% Toolbox dir
%------------------------------------------------------------
disp('Toolbox');
TbxSrcDir = [SourceDir filesep 'easyspin'];

% Asserting that c files do not contain // comments (which are not
% supported by all C compilers.)
fprintf('  Checking *.c files for absence of // comments.');
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
fprintf('ok\n');

disp('  Generating toolbox folder.');
mkdir(BuildDir,'easyspin');

disp('  Copying toolbox files');
copyfile(TbxSrcDir,TbxDir);

%{
% Remove Mercurial repository and files
disp('  Removing source control files.');
if exist([TbxDir filesep '.hg'],'dir')
  rmdir([TbxDir filesep '.hg'],'s');
end
if exist([TbxDir filesep '.hgignore'],'file')
  delete([TbxDir filesep '.hgignore']);
end
if exist([TbxDir filesep '.hgtags'],'file')
  delete([TbxDir filesep '.hgtags']);
end
%}

% Release tag and expiry date
disp('  Setting release information and expiry date.');
replacestr([TbxDir filesep 'info.xml'],'$ReleaseID$',ReleaseID);
replacestr([TbxDir filesep 'easyspininfo.m'],'$ReleaseID$',ReleaseID);
replacestr([TbxDir filesep 'easyspininfo.m'],'$ReleaseDate$',ReleaseDate);
replacestr([TbxDir filesep 'easyspininfo.m'],'$ExpiryDate$',ExpiryDate);
replacestr([TbxDir filesep 'eschecker.m'],'888888',num2str(datenum(ExpiryDate)));
replacestr([TbxDir filesep 'eschecker.m'],'999999',num2str(datenum(HorizonDate)));

%---------------------------------------------------------------------
% P-coding
%---------------------------------------------------------------------
% Copy everything to a third directory for pcoding
% This makes sure .m files are not newer than .p files
%  (otherwise Matlab complains about potentially obsolete p-files)
disp('  Copying to a new folder for p-coding.');
mkdir(BuildDir,'easyspinpcode');
copyfile(TbxDir,TbxPcodeDir);

disp('  Extracting help from all m-files.');
extracthelp(TbxDir);

disp('  Deleting m-files in private subdirectory.');
delete([TbxDir filesep 'private' filesep '*.m']);

curdir = pwd;
disp('  P-coding all m-files.');
cd(TbxDir);
pcode([TbxPcodeDir filesep '*.m']);
cd([TbxDir filesep 'private']);
pcode([TbxPcodeDir filesep 'private' filesep '*.m']);
cd(curdir);

rmdir(TbxPcodeDir,'s');

%------------------------------------------------------------
% Examples
%------------------------------------------------------------
disp('Examples');

disp('  generating examples dir and copying files');
mkdir(BuildDir,'examples');
copyfile([SourceDir filesep 'examples'],ExmplDir);

disp('  generating examples html file');
perl('mkexamples.pl');


%------------------------------------------------------------
% Documentation
%------------------------------------------------------------
disp('Documentation');
disp('  generating documentation folder');
mkdir(BuildDir,'documentation');

DocSrcDir = [SourceDir filesep 'docs'];
disp('  copying files');
copyfile(DocSrcDir,DocDir);

% Replace all $ReleaseID$ in html files with actual release ID
disp('  inserting release ID into documentation');
docFiles = dir([DocDir filesep '*.html']);
for iFile = 1:numel(docFiles)
  replacestr([DocDir filesep docFiles(iFile).name],'$ReleaseID$',ReleaseID);
end

%{
% Remove Mercurial repository and files
disp('  removing source control files');
if exist([DocDir filesep '.hg'],'dir')
  rmdir([DocDir filesep '.hg'],'s');
end
if exist([DocDir filesep '.hgignore'],'file')
  delete([DocDir filesep '.hgignore']);
end
if exist([DocDir filesep '.hgtags'],'file')
  delete([DocDir filesep '.hgtags']);
end
%}

%------------------------------------------------------------
% Packaging
%------------------------------------------------------------
disp('Packaging');
% package for public release
if ZipDestDir(end)==filesep, ZipDestDir(end) = []; end
zipFile = [ZipDestDir filesep 'easyspin-' BuildID '.zip'];
zipFileShort = [ZipDestDir filesep 'easyspin-' ReleaseID '.zip'];

disp(['  packing zip file ' zipFile '.']);
if exist(zipFile,'file'),
  try
    delete(zipFile);
  catch
    error('Cannot delete zip file.');
  end
end
zip(zipFile,BuildDir);

if ~betaVersion
  disp(['  copying to ' zipFileShort '.']);
  copyfile(zipFile,zipFileShort,'f');
end

% entire development tree
if packDevtree
    zipFile = [ZipDestDir filesep 'easyspin-' BuildID '-devel.zip'];
    disp(['  packing zip file ' zipFile ' .']);
    if exist(zipFile,'file'),
        try
            delete(zipFile);
        catch
            error('Cannot delete zip file.');
        end
    end
    DevelopmentDir = SourceDir;
    zip(zipFile,DevelopmentDir);
end

%------------------------------------------------------------
% Clean-up
%------------------------------------------------------------
disp('Removing build tree.');
rmdir(BuildDir,'s');

save lastbuildid BuildNumber

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
Lines = textscan(fid,'%s','whitespace','','delimiter','\n','bufsize',100000);
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
