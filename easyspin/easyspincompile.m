function easyspincompile

warning('off','MATLAB:oldPfileVersion');

disp('EasySpin compilation');

% Determine directory containing mex source files
%-----------------------------------------------------
esPath = fileparts(which(mfilename));
mexDirectory = [esPath filesep 'private'];
disp(['  directory: ' mexDirectory]);
disp(['  version: ' version]);
% Change directory
olddir = pwd;
cd(mexDirectory);


% List of *.c files
%-----------------------------------------------------
SourceFiles = dir('*.c');



% Determine architecture and set mex options
%-----------------------------------------------------
if strcmp(mexext,'mexa64') || strcmp(mexext,'mexw64') || strcmp(mexext,'mexmaci64');
  nBits = 64;
  mexoptions = '-Dbit64';
else
  nBits = 32;
  mexoptions = '-Dbit32';
end
fprintf('  mex extension: %s, %d-bit\n',mexext,nBits);


% Compile and link mex files
%-----------------------------------------------------
fprintf('  compiling..');
for f = 1:numel(SourceFiles)
  fprintf('.');
  mex(mexoptions,SourceFiles(f).name);
end
disp('done.');
% Restore old directory
cd(olddir);

% A hack that was needed at the EasySpin workshop at Cornell 2007,
% no idea why. Needed when the C compiler for mex files was not
% configured.
clear functions
