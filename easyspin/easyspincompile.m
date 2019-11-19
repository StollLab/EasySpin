function easyspincompile

warning('off','MATLAB:oldPfileVersion');

disp('EasySpin compilation');

% Determine directory containing mex source files
%-------------------------------------------------------------------------------
esPath = fileparts(which(mfilename));
mexDirectory = [esPath filesep 'private'];
disp(['  directory: ' mexDirectory]);
disp(['  version: ' version]);

% Change directory
olddir = pwd;
cd(mexDirectory);


% Determine architecture and set mex options
%-------------------------------------------------------------------------------
if strcmp(mexext,'mexa64') || strcmp(mexext,'mexw64') || strcmp(mexext,'mexmaci64')
  nBits = 64;
  mexoptions = '-Dbit64';
else
  nBits = 32;
  mexoptions = '-Dbit32';
end
mexoptions = {mexoptions,'-silent'};
fprintf('  mex extension: %s, %d-bit\n',mexext,nBits);


% Determine mex configuration
%-------------------------------------------------------------------------------
cc = mex.getCompilerConfigurations('C','Installed');
if numel(cc)==0
  error('MEX is not configured for C. Run mex -setup C.');
else
  fprintf('  mex C compiler: %s\n',cc(1).Name); 
end


% Get list of *.c files
%-------------------------------------------------------------------------------
SourceFiles = dir('*.c');
nFiles = numel(SourceFiles);


% Compile and link mex files
%-------------------------------------------------------------------------------
fprintf('  compiling %d c-files...\n',nFiles);
for f = 1:nFiles
  fprintf('  (%d/%d) %-25s ',f,nFiles,SourceFiles(f).name);
  try
    mex(SourceFiles(f).name,mexoptions{:});
    fprintf('  complete\n');
    ok(f) = true;
  catch
    fprintf('  failed\n');
    disp(lasterr);
    ok(f) = false;
  end
end

cd(olddir);

% A hack that was needed at the EasySpin workshop at Cornell 2007, no idea why.
% it works. Needed when the C compiler for mex files was not configured.
clear functions
