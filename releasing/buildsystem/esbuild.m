%===============================================================================
% EasySpin build script
%===============================================================================
function esbuild()

fprintf('Building EasySpin.\n');

% Configuration
%-------------------------------------------------------------------------------
if exist('esbuild_config.m','file')
  fprintf('Reading configuration file...');
  try
    esbuild_config
  catch
    error('Error executing esbuild_config.m.');
  end
  fprintf('ok\n');
else
  error('esbuild_config.m not found.');
end

fprintf('Configuration:\n')
fprintf('  releaseID       %s\n',releaseID);
fprintf('  releaseChannel  %s\n',releaseChannel);
fprintf('  monthsToExpiry  %d\n',monthsToExpiry);
fprintf('  sourceDir       %s\n',sourceDir);
fprintf('  zipDestDir      %s\n',zipDestDir);

if releaseID(1)=='v'
  error('Release IDs should not start with v.');
end

% Set release and expiry date
nowDate = datetime('now','Format','yyyy-MM-dd');
releaseDate = char(nowDate);
expiryDate = nowDate + calendarDuration(0,monthsToExpiry,0);
expiryDate = char(expiryDate);

fprintf('Current directory:    %s\n',pwd);
fprintf('Source directory: %s\n',sourceDir);

fprintf('Checking for source subdirectories...');

if ~exist([sourceDir filesep 'easyspin'],'dir')
  error('Could not find EasySpin source directory easyspin/.');
end
if ~exist([sourceDir filesep 'examples'],'dir')
  error('Could not find EasySpin source directory examples/.');
end

fprintf(' ok\n');

fprintf('Creating temporary build directory...');

buildDir = [tempdir,'easyspin-',releaseID];

if exist(buildDir,'dir')
  ok = rmdir(buildDir,'s');
  if ~ok
    error('Could not remove existing temporary build directory %s.',buildDir);
  end
end
if exist(buildDir,'file')
  delete(buildDir);
end

ok = mkdir(buildDir);
if ~ok
  error('Could not create temporary build directory %s.',buildDir);
end

fprintf(' ok\n');
fprintf('  Temp build directory: %s\n',buildDir);

tbxDir = [buildDir filesep 'easyspin'];
exmplDir = [buildDir filesep 'examples'];
docDir = [buildDir filesep 'documentation'];

tbxPcodeDir = [buildDir filesep 'easyspinpcode'];

% Toolbox folder
%-------------------------------------------------------------------------------
disp('Toolbox');
tbxSrcDir = [sourceDir filesep 'easyspin'];

% Assert that .c files do not contain // comments (which are not
% supported by strict ANSI-C compilers, but are supported by the C99
% standard .)
enforceANSICcomments = false;
if enforceANSICcomments
  fprintf('  Checking *.c files for absence of // comments...');
  cfiles = dir([tbxSrcDir filesep 'private' filesep '*.c']);
  for k = 1:numel(cfiles)
    Lines = textread([tbxSrcDir filesep 'private' filesep cfiles(k).name],'%s','whitespace','\n');
    if any(contains(Lines,'//'))
      error('Found // comment in file %s',['private' filesep cfiles(k).name]);
    end
  end
  fprintf(' ok\n');
end

fprintf('  Generating toolbox directory...');
mkdir(buildDir,'easyspin');
fprintf(' ok\n');

fprintf('  Copying toolbox files... ');
copyfile(tbxSrcDir,tbxDir);
fprintf(' ok\n');


% Release information
%-------------------------------------------------------------------------------
fprintf('  Setting release information... ');

esinfoFile = [tbxDir filesep 'private' filesep 'easyspin_info.m'];
replacestr(esinfoFile,'$ReleaseID$',releaseID);
replacestr(esinfoFile,'$ReleaseChannel$',releaseChannel);
replacestr(esinfoFile,'$ReleaseDate$',releaseDate);
replacestr(esinfoFile,'$ExpiryDate$',expiryDate);
fprintf(' ok\n');


% P-coding
%-------------------------------------------------------------------------------
% Copy everything to a temp directory for pcoding
% This makes sure .m files are not newer than .p files
%  (otherwise Matlab complains about potentially obsolete p-files)
fprintf('  Copying to a new directory for p-coding...');
mkdir(buildDir,'easyspinpcode');
copyfile(tbxDir,tbxPcodeDir);
fprintf(' ok\n');

fprintf('  Extracting help from all m-files...');
extracthelp(tbxDir);
fprintf(' ok\n');

fprintf('  Deleting m-files in private subdirectory...');
delete([tbxDir filesep 'private' filesep '*.m']);
fprintf(' ok\n');

% Pcode format changed several times: R2007b and R2022a
% We use Pcode format version R2007b-R2021b
fprintf('  P-coding all m-files...');
MatlabVersion = sscanf(version,'%d.%d',2);
addPcodeVersionFlag = MatlabVersion(1)>9 || (MatlabVersion(1)==0 && MatlabVersion(2)>=12);  % 9.12 = R2022a
curdir = pwd;
cd(tbxDir);
if addPcodeVersionFlag
  pcode([tbxPcodeDir filesep '*.m'],'-R2007b');
else
  pcode([tbxPcodeDir filesep '*.m']);
end
cd([tbxDir filesep 'private']);
if addPcodeVersionFlag
  pcode([tbxPcodeDir filesep 'private' filesep '*.m'],'-R2007b');
else
  pcode([tbxPcodeDir filesep 'private' filesep '*.m']);
end
cd(curdir);
rmdir(tbxPcodeDir,'s');
fprintf(' ok\n');



% Examples
%-------------------------------------------------------------------------------
disp('Examples');

fprintf('  creating examples dir and copying files...');
mkdir(buildDir,'examples');
copyfile([sourceDir filesep 'examples'],exmplDir);
fprintf(' ok\n');


% Documentation
%-------------------------------------------------------------------------------
disp('Documentation');

fprintf('  creating documentation folder...');
mkdir(buildDir,'documentation');
fprintf(' ok\n');

fprintf('  copying files...');
docSrcDir = [sourceDir filesep 'documentation'];
copyfile(docSrcDir,docDir);
fprintf(' ok\n');

% Replace all $ReleaseID$ in HTML files with actual release ID
fprintf('  inserting release ID into documentation...');
docFiles = dir([docDir filesep '*.html']);
for iFile = 1:numel(docFiles)
  replacestr([docDir filesep docFiles(iFile).name],'$ReleaseID$',releaseID);
end
fprintf(' ok\n');


% Packaging
%-------------------------------------------------------------------------------
disp('Packaging');
if zipDestDir(end)==filesep, zipDestDir(end) = []; end
zipFileName = [zipDestDir filesep 'easyspin-' releaseID '.zip'];

fprintf('  packing zip file %s...',zipFileName);
zip(zipFileName,buildDir);
fprintf(' ok\n');


% Clean-up
%-------------------------------------------------------------------------------
fprintf('Removing temporary build folder...');
rmdir(buildDir,'s');
fprintf(' ok\n');

disp('Done!');

end


%===============================================================================

% extracthelp  Extracts help from m files
%
%    extracthelp(pDir)
%
%    Extract the help from all m-files in the directory
%    pDir, i.e. removes all the code.

function extracthelp(pDir)

% Get list of all m-files in directory
mFiles = dir([pDir filesep '*.m']);
mFiles = sort({mFiles.name});
for k = 1:length(mFiles)
  mFiles{k} = [pDir filesep mFiles{k}];
end

% Copy help lines from dir/file.m to dir2/file.m
for k = 1:length(mFiles)
  infile = fopen(mFiles{k},'r');
  commentTag = '%';
  ok = true;
  lineStr = {};
  while ok
    thisLine = fgetl(infile);
    if isempty(thisLine) || all(thisLine==-1); break; end
    ok = thisLine(1)==commentTag;
    if ok, lineStr{end+1} = thisLine; end
  end
  fclose(infile);

  outfile = fopen(mFiles{k},'w');
  for l = 1:numel(lineStr)
    fprintf(outfile,'%s\n',lineStr{l});
  end
  fclose(outfile);
end

end

%===============================================================================
% replacestr  Replace all occurances of a string in a file
function replacestr(fname,oldString,newString)
txt = fileread(fname);
newtxt = replace(txt,oldString,newString);
f = fopen(fname,'w');
fwrite(f,newtxt);
fclose(f);
end
