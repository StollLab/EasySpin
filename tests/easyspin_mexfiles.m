function err = test()

% Check for presence of all mex files for all platforms
%===============================================================================

olddir = pwd;
cd('../easyspin/private/');

SourceFiles = dir('*.c');

% Assume that all *.c files give *.mex* files
mexFiles = SourceFiles;
nFiles = numel(mexFiles);

ext = {'.mexw64','.mexmaci64','.mexa64'};
errormsg = '';
nMissing = 0;
for k = 1:nFiles
  fn = mexFiles(k).name(1:end-2);
  for e = 1:3
    fnmex = [fn ext{e}];
    if ~exist(fnmex,'file')
      nMissing = nMissing+1;
      errormsg = [errormsg sprintf('    %s\n',fnmex)];
    end
  end
end
cd(olddir);

err = nMissing~=0;

if err
  errormsg = [sprintf('  %d mex files are missing:\n',nMissing) errormsg];
  error(errormsg);
end
