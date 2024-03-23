% estest    Unit testing engine for EasySpin
%
%   Usage:
%     estest            run all tests
%     estest asdf       run all tests whose name starts with asdf
%     estest asdf d     run all tests whose name starts with asdf and
%                       display results
%     estest asdf r     evaluate all tests whose name starts with asdf and
%                       recalculate and store regression data
%     estest asdf t     evaluate all tests whose name starts with asdf and
%                       report timings
%     estest asdf l     evaluate all tests whose name starts with asdf and
%                       report all lines of code not covered by the tests
%
%   Either the command syntax as above or the function syntax, e.g.
%   estest('asdf','t'), can be used.
%
%   Run all tests including timings:   estest('*_','t')
%
%   All test files must have an underscore _ in their filename.

function out = estest(TestName,params)

% Check whether EasySpin is on the MATLAB path
EasySpinPath = fileparts(which('easyspin'));
if isempty(EasySpinPath)
  error('EasySpin is not on the MATLAB path!');
end

fid = 1; % output to command window

if nargin<1
  TestName = '';
end
if nargin<2
  params = '';
end

Opt.Display = any(params=='d');
Opt.Regenerate = any(params=='r');
Opt.Verbosity = Opt.Display;

displayTimings = any(params=='t');
runCodeCoverageAnalysis = any(params=='l');

if Opt.Display && displayTimings
  error('Cannot plot test results and report timings at the same time.');
end

if any(TestName=='_')
  FileMask = [TestName '*.m'];
else
  FileMask = [TestName '*_*.m'];
end

FileList = dir(FileMask);

if numel(FileList)==0
  fprintf('No test functions matching the pattern %s\n',FileMask);
  return
end

TestFileNames = sort({FileList.name});

fprintf(fid,'=======================================================================\n');
fprintf(fid,'EasySpin test set                      %s\n(MATLAB %s)\n',char(datetime),version);
fprintf(fid,'EasySpin location: %s\n',EasySpinPath);
fprintf(fid,'=======================================================================\n');
fprintf(fid,'Display: %d, Regenerate: %d, Verbosity: %d\n',...
  Opt.Display,Opt.Regenerate, Opt.Verbosity);
fprintf(fid,'-----------------------------------------------------------------------\n');

% test outcome codes:
%    0   test passed
%   +1   test failed
%   +2   test crashed
%   +3   not tested

OutcomeStrings = {'pass','failed','crashed','not tested'};

% Get path to Easyspin functions folder
path = fileparts(which('estest'));
path = path(1:end-length('\tests'));

% List all the functions, including private
Files = dir(fullfile(path,'easyspin','*.m'));
executedLines = repmat({[]},length(Files),1);

for iTest = 1:numel(TestFileNames)
  
  thisTestName = TestFileNames{iTest}(1:end-2);

  if Opt.Display
    clf
    set(gcf,'Name',thisTestName);
    drawnow
  end
  
  % Load, or regenerate, reference data
  refdata = [];
  TestDataFile = ['data/' thisTestName '.mat'];
  if exist(TestDataFile,'file')
    if Opt.Regenerate
      delete(TestDataFile);
      refdata = [];
    else
      try
        refdata = load(TestDataFile,'data');
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
  time_used(iTest) = toc;

  if Opt.Display
    if iTest<numel(TestFileNames), pause; end
  end  
  
  % Retrieve profiler summary and turn profiler off
  if runCodeCoverageAnalysis
    p = profile('info');
    profile off
    
    % Make list of all profiled function calls
    executedFcns = {p.FunctionTable(:).CompleteName};
    % Analyze code coverage of each API function
    for n = 1:length(Files)
      FcnName = Files(n).name;
      pos = find(contains(executedFcns,FcnName));
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
    save(TestDataFile,'data');
  end
  
  testResults(iTest).outcome = testOutcome;
  testResults(iTest).name = thisTestName;
  testResults(iTest).errorData = errorInfo;
  
  outcomeStr = OutcomeStrings{testResults(iTest).outcome+1};
  
  if ~all(ok(:)) && ~displayTimings
    outcomeStr = [outcomeStr '  ' num2str(find(~ok(:).'))];
  end
  
  if ~isempty(data)
    typeStr = 'regression';
  else
    typeStr = 'direct';
  end
  
  if displayTimings
    timeStr = sprintf('%0.3f seconds',time_used(iTest));
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
  fprintf(fid,'Code Coverage Analysis \n');
  fprintf(fid,'-----------------------------------------------------------------------\n');
  
  TotalCovered = 0;
  TotalRunnable = 0;
  % Analyze code coverage of each function
  for n = 1:length(Files)
    FcnName = Files(n).name;
    Path = Files(n).folder;
    RunnableLines = callstats('file_lines',fullfile(Path,FcnName));
    if isempty(RunnableLines)
      RunnableLines = 0;
    end
    Executed = unique(executedLines{n});
    Covered = length(Executed);
    TotalCovered = TotalCovered + Covered;
    Runnable = length(unique(RunnableLines));
    Code = fileread(FcnName);
    % Account for lines missed by callstats
    if Covered > Runnable
        TotalRunnable = TotalRunnable + Covered;
    else
        TotalRunnable = TotalRunnable + Runnable;
    end
    if params =='l'
      Missed = RunnableLines;
      for k=1:length(Executed)
        Missed(RunnableLines==Executed(k)) = NaN;
      end
      Missed(isnan(Missed)) = [];
    end
    % Account for unreachable end-statements or lines
    MissedEnds = length(strfind(Code,'error')) + length(strfind(Code,'return')) ...
      + length(strfind(Code,'break')) + length(strfind(Code,'fprintf')) + length(strfind(Code,'plot'));
    if Runnable - Covered <= MissedEnds
      Covered = Runnable;
      Missed = [];
    end
    Coverage = 100*Covered/Runnable;
    % Print to command window
    if (~isempty(TestName) && Coverage~=0) || isempty(TestName)
        fprintf('%-20s%-18s%5.1f%% %18s  %s\n',FcnName,' ',Coverage,'Lines missing:',mat2str(Missed))
    end
  end
  TotalCoverage = TotalCovered/TotalRunnable*100;
  if isempty(TestName)
    fprintf('Total code coverage: %3.2f%%\n',TotalCoverage)
  end
end

allOutcomes = [testResults.outcome];

% Display timings of slowest tests
if displayTimings
  fprintf(fid,'-----------------------------------------------------------------------\n');
  fprintf(fid,'Total test time:                        %7.3f seconds\n',sum(time_used));
  fprintf(fid,'Slowest tests:\n');
  [time,iTest] = sort(time_used,'descend');
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

% Return output if desired
if nargout==1
  out.Results = testResults;
  out.outcomes = allOutcomes;
end

return
