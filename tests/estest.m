% estest    Unit test runner for EasySpin
%
%   Usage:
%     estest            run all tests
%     estest asdf       run all tests whose names start with asdf
%     estest asdf d     run all tests whose names start with asdf and
%                       display results
%     estest asdf r     run all tests whose names start with asdf and
%                       recalculate and store regression data
%     estest asdf t     run all tests whose names start with asdf and
%                       report timings
%     estest asdf c     run all tests whose name starts with asdf and
%                       report all lines of code not covered by the tests
%
%   Either the command syntax, as above, or the function syntax, e.g.
%   estest('asdf','t'), can be used.
%
%   Run all tests including timings:   estest('','t')
%
%   To be recognized as a test function, the file name must contain an
%   underscore (_) in their filename. Files without _ are ignored.
%
%   A test function must return one output, ok, which can be a scalar or
%   an array of true/false. True indicates that the test has passed. An
%   array indicates a series of subtests.
%
%   Optionally, a structure is passed to the test function that contains
%   the following fields:
%
%   Opt.Display     true/false - whether the test function should plot/print
%   Opt.Regenerate  true/false - whether the test function should
%                    regenerate reference data
%   Opt.Verbosity   true/false - is passed on to EasySpin functions for
%                    additional logging

function out = estest(testName,params)

% Check whether EasySpin is on the MATLAB path
EasySpinPath = fileparts(which('easyspin'));
if isempty(EasySpinPath)
  error('EasySpin is not on the MATLAB path!');
end

fid = 1;  % output to command window

if nargin<1
  testName = '';
end

if nargin<2
  params = '';
end

% Options to pass along to test functions
Opt.Display = any(params=='d');
Opt.Regenerate = any(params=='r');
Opt.Verbosity = double(Opt.Display);

displayTimings = any(params=='t');
runCodeCoverageAnalysis = any(params=='c');

if Opt.Display && displayTimings
  error('Cannot plot test results and report timings at the same time.');
end

% Assemble list of file names of tests to be run
if ischar(testName)
  if any(testName=='*')
    error('Dont''t use * in first input argument.')
  end
  % Get list of m files in folder, remove all file without _
  fileMask = [testName '*.m'];
  fileList = dir(fileMask);
  for i = numel(fileList):-1:1
    if ~any(fileList(i).name=='_')
      fileList(i) = [];
    end
  end
  if numel(fileList)==0
    fprintf('No test functions matching the pattern %s\n',fileMask);
    return
  end
  testFileNames = {fileList.name}.';
elseif iscell(testName)
  if numel(testName)==0
    fprintf('No tests to run.\n');
    return
  end
  testName = cellstr(testName);  % convert strings to char arrays (needed by sort)
  testFileNames = cellfun(@(x)[x '.m'],testName,'UniformOutput',false);
else
  error('First input must be a string/character array, or a cell array of such.');
end
testFileNames = sort(testFileNames);

fprintf(fid,'=======================================================================\n');
fprintf(fid,'EasySpin test set                      %s\n(MATLAB %s)\n',char(datetime),version);
fprintf(fid,'EasySpin folder: %s\n',EasySpinPath);
fprintf(fid,'=======================================================================\n');
fprintf(fid,'Display: %d, Regenerate: %d, Verbosity: %d\n',...
  Opt.Display,Opt.Regenerate,Opt.Verbosity);
fprintf(fid,'-----------------------------------------------------------------------\n');

% test outcome codes:
%    0   test passed
%   +1   test failed
%   +2   test crashed
%   +3   not tested

outcomeStrings = {'pass','failed','crashed','not tested'};

% List all the functions, including private
if runCodeCoverageAnalysis
  % Get path to Easyspin functions folder
  path = fileparts(which('estest'));
  path = path(1:end-length('\tests'));
  Files = dir(fullfile(path,'easyspin','*.m'));
  executedLines = repmat({[]},length(Files),1);
end

nTests = numel(testFileNames);
timeElapsed = zeros(1,nTests);
testResults(nTests) = struct;
for iTest = 1:nTests
  
  thisTestName = testFileNames{iTest}(1:end-2);

  if Opt.Display
    clf
    set(gcf,'Name',thisTestName);
    drawnow
  end
  
  % Load, or regenerate, reference data
  refdata = [];
  testDataFile = ['data/' thisTestName '.mat'];
  if exist(testDataFile,'file')
    if Opt.Regenerate
      delete(testDataFile);
      refdata = [];
    else
      try
        refdata = load(testDataFile,'data');
        refdata = refdata.data;
      catch
        error('Could not load data for test ''%s''.',thisTestName);
      end
    end
  end
  
  % Clear and start profiler
  if runCodeCoverageAnalysis
    profile clear
    profile on
  end
  
  % Run test, catch any errors
  testFcn = str2func(thisTestName);
  nArgsOut = nargout(testFcn);
  nArgsIn = nargin(testFcn);
  usesStoredData = nArgsIn==2 && nArgsOut==2;
  tic
  try
    if usesStoredData
      if nArgsIn<2, error('2 inputs are needed.'); end
      [ok,data] = testFcn(Opt,refdata);
    else
      if nArgsIn==0
        ok = testFcn();
        data = [];
      else
        ok = testFcn(Opt);
        data = [];
      end
    end
    if isempty(ok)
      testOutcome = 3; % not tested
    else
      if all(ok)
        testOutcome = 0; % test passed
      else
        testOutcome = 1; % test failed
      end
    end
    errorInfo = [];
    errorStr = '';
  catch exception
    ok = false;
    data = [];
    testOutcome = 2; % test crashed
    errorInfo = exception;
    errorStr = getReport(errorInfo);
    errorStr = ['    ' regexprep(errorStr,'\n','\n    ') newline];
  end
  timeElapsed(iTest) = toc;

  if Opt.Display
    if iTest<numel(testFileNames)
      pause;
    end
  end  
  
  % Retrieve profiler summary and turn profiler off
  if runCodeCoverageAnalysis
    p = profile('info');
    profile off
    
    % Make list of all profiled function calls
    executedFcns = {p.FunctionTable(:).CompleteName};
    % Analyze code coverage of each API function
    for n = 1:length(Files)
      fcnName = Files(n).name;
      pos = find(contains(executedFcns,fcnName));
      if ~isempty(pos)
        % initialize containers
        for i = 1:length(pos)
          % get executed lines in profiler
          tmp = p.FunctionTable(pos(i)).ExecutedLines;
          container = executedLines{n};
          container(end+1:end+length(tmp(:, 1))) = tmp(:, 1);
          executedLines{n} = container;
        end
      end
    end
  end
  
  saveTestData = usesStoredData && Opt.Regenerate;
  if saveTestData
    save(testDataFile,'data');
  end
  
  testResults(iTest).outcome = testOutcome;
  testResults(iTest).name = thisTestName;
  testResults(iTest).errorData = errorInfo;
  
  outcomeStr = outcomeStrings{testResults(iTest).outcome+1};
  if ~all(ok(:)) && ~displayTimings
    outcomeStr = [outcomeStr '  ' num2str(find(~ok(:).'))];
  end
  
  if ~isempty(data)
    typeStr = 'regression';
  else
    typeStr = 'direct';
  end
  
  if displayTimings
    timeStr = sprintf('%0.3f seconds',timeElapsed(iTest));
  else
    timeStr = '';
  end

  nameStr = testResults(iTest).name;
  str = sprintf('%-47s  %-12s%-8s%s\n%s',...
       nameStr,typeStr,outcomeStr,timeStr,errorStr);
  str(str=='\') = '/';
  
  nBlanks = max(47-length(nameStr),0);
  nameStrLink = sprintf('<a href="matlab: edit %s">%s</a>%s',nameStr,nameStr,repmat(' ',1,nBlanks));
  strLink = sprintf('%s  %-12s%-8s%s\n%s',...
       nameStrLink,typeStr,outcomeStr,timeStr,errorStr);
  strLink(strLink=='\') = '/';
  
  testResults(iTest).msg = str;
  testResults(iTest).msgLink = strLink;
  
  fprintf(fid,strLink);
end

% Display results of code coverage analysis
if runCodeCoverageAnalysis
  
  fprintf(fid,'-----------------------------------------------------------------------\n');
  fprintf(fid,'Code coverage analysis \n');
  fprintf(fid,'-----------------------------------------------------------------------\n');
  
  totalCovered = 0;
  totalRunnable = 0;
  % Analyze code coverage of each function
  for n = 1:length(Files)
    fcnName = Files(n).name;
    Path = Files(n).folder;
    runnableLines = callstats('file_lines',fullfile(Path,fcnName));
    if isempty(runnableLines)
      runnableLines = 0;
    end
    executed = unique(executedLines{n});
    covered = length(executed);
    totalCovered = totalCovered + covered;
    runnable = length(unique(runnableLines));
    Code = fileread(fcnName);
    % Account for lines missed by callstats
    if covered > runnable
        totalRunnable = totalRunnable + covered;
    else
        totalRunnable = totalRunnable + runnable;
    end
    if params =='l'
      missed = runnableLines;
      for k=1:length(executed)
        missed(runnableLines==executed(k)) = NaN;
      end
      missed(isnan(missed)) = [];
    end
    % Account for unreachable end-statements or lines
    missedEnds = length(strfind(Code,'error')) + length(strfind(Code,'return')) ...
      + length(strfind(Code,'break')) + length(strfind(Code,'fprintf')) + length(strfind(Code,'plot'));
    if runnable - covered <= missedEnds
      covered = runnable;
      missed = [];
    end
    coverage = 100*covered/runnable;
    % Print to command window
    if (~isempty(testName) && coverage~=0) || isempty(testName)
        fprintf('%-20s%-18s%5.1f%% %18s  %s\n',fcnName,' ',coverage,'Lines missing:',mat2str(missed))
    end
  end
  totalCoverage = totalCovered/totalRunnable*100;
  if isempty(testName)
    fprintf('Total code coverage: %3.2f%%\n',totalCoverage);
  end
end

allOutcomes = [testResults.outcome];

% Display timings of slowest tests
if displayTimings
  fprintf(fid,'-----------------------------------------------------------------------\n');
  fprintf(fid,'Total test time:                        %7.3f seconds\n',sum(timeElapsed));
  fprintf(fid,'Slowest tests:\n');
  [time,iTest] = sort(timeElapsed,'descend');
  for q = 1:min(10,numel(time))
    fprintf(fid,'%-36s    %7.3f seconds\n',testResults(iTest(q)).name,time(q));
  end
end

% Display all tests that failed or crashed
if any(allOutcomes==1) || any(allOutcomes==2)
  fprintf(fid,'-----------------------------------------------------------------------\n');
  for iTest = find(allOutcomes)
    fprintf(fid,testResults(iTest).msgLink);
  end
end

nPasses = sum(allOutcomes==0);
nFailures = sum(allOutcomes==1);
nCrashes = sum(allOutcomes==2);
fprintf(fid,'-----------------------------------------------------------------------\n');
msg = sprintf('%d passes, %d failures, %d crashes\n',nPasses,nFailures,nCrashes);
fprintf(fid,msg);
fprintf(fid,'-----------------------------------------------------------------------\n');

% Assemble list of all failed and crashed tests and provide link to rerun them.
if nFailures+nCrashes>0
  testList = cell(1,nFailures+nCrashes);
  idx = 1;
  for iTest = find(allOutcomes)
    testList{idx} = ['''' testResults(iTest).name ''''];
    idx = idx + 1;
  end
  testList = ['{' strjoin(testList,',') '}'];
  fprintf(fid,'<a href="matlab: estest(%s)">rerun failed and crashed tests</a>\n',testList);
end

% Return output if desired
if nargout==1
  out.Results = testResults;
  out.outcomes = allOutcomes;
end

end
