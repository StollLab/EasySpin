% estestexamples    Check EasySpin examples for crashes
%
%   Usage:
%     estestexamples              % run all examples
%     estestexamples foldername   % run examples in given folder
%

function out = estestexamples(examplefolder)

% Check whether EasySpin is on the MATLAB path
EasySpinPath = fileparts(which('easyspin'));
if isempty(EasySpinPath)
  error('EasySpin is not on the MATLAB path!');
end

% Get path to examples
ExamplesPath = [EasySpinPath(1:end-8) 'examples'];

if nargin<1
  ExamplesFolders = dir(ExamplesPath);
  ExamplesFolders(strcmp({ExamplesFolders.name},'.')) = [];
  ExamplesFolders(strcmp({ExamplesFolders.name},'..')) = [];
else
  if ~isfolder(fullfile(ExamplesPath,examplefolder))
    error('Provide a valid example folder as input or run without arguments to test all examples.')
  end
  ExamplesFolders(1).name = examplefolder;  
end

% Loop over folders, find and run all example *.m files
% ------------------------------------------------------------------------------------------------
s = 0;
for i = 1:numel(ExamplesFolders)

  f = dir(fullfile(ExamplesPath,ExamplesFolders(i).name,'*.m'));

  for j = 1:numel(f)

    % Run example and check if it runs without errors
    clear exampleend errormsg
    examplestart = datetime;
    save('workspacetmp.mat') % save workspace to recover variables if example contains clear
    try
      run(fullfile(f(j).folder,f(j).name))
      errormsg = [];
    catch exception
      errormsg = exception.message;
    end
    exampleend = datetime;
    clearvars -except errormsg exampleend % reset workspace
    load('workspacetmp.mat') % reload variables

    % Save result in example structure
    s = s+1;
    example(s).folder = ExamplesFolders(i).name;
    example(s).name = f(j).name;
    example(s).duration = seconds(exampleend-examplestart);
    if isempty(errormsg)
      example(s).error = 0;
      example(s).errormsg = '';
    else
      example(s).error = 1;
      example(s).errormsg = errormsg;
    end

  end

end
delete('workspacetmp.mat')

% Summary
% ------------------------------------------------------------------------------------------------
% List examples and duration
clc
fprintf('\n')
fprintf('List of examples and duration (sorted from slowest to fastest):\n')
fprintf('-------------------------------------------------------------------------\n')
[~,ind] = sort([example.duration],'descend');
for i = ind
  if example(i).error
  fprintf('%-52s: %1.3f s Error: %s \n',fullfile(example(i).folder,example(i).name),example(i).duration,example(i).errormsg)
  else
  fprintf('%-52s: %1.3f s \n',fullfile(example(i).folder,example(i).name),example(i).duration)
  end
end

% Check if any examples crashed
errorind = find([example.error]);
ncrashes = numel(errorind);
fprintf('\n')
fprintf('Summary\n')
fprintf('-------------------------------------------------------------------------\n')
fprintf('Total example test time: %1.1f min\n',sum([example.duration])/60)
fprintf('%i passes, %i crashes\n',numel(example)-ncrashes,ncrashes)
fprintf('-------------------------------------------------------------------------\n')

% List crashed examples and error messages
if ~isempty(errorind)
  fprintf('\n')
  fprintf('Examples with error messages:\n')
  fprintf('-------------------------------------------------------------------------\n')
  for i = errorind
    fprintf('%-52s: %s \n',fullfile(example(i).folder,example(i).name),example(i).errormsg)
  end
end

% Return output if desired
if nargout==1
  out = example;
end