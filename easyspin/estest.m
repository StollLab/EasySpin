% estest    Unit test runner for EasySpin
%
%   Usage:
%     estest                  run all tests
%     estest *                run all tests
%     estest pepper_c2h       run only the test pepper_c2h
%     estest pepper*          run all tests whose names start with pepper
%     estest *crystal*        run all tests whose names contain crystal
%     estest pepper* garlic*  run all tests matching any of the patterns
%     estest -t               run all tests and report timings
%     estest pepper* -t       run all tests starting with pepper and report
%                             timings
%
%   Arguments starting with - are options, all others are test name
%   patterns. A pattern matches exactly, unless it contains the wildcard *.
%   Options can be given in any position and can be combined (e.g. -tc).
%
%   Do not use a lone * followed by further arguments in command syntax:
%   MATLAB parses "estest * -t" as the expression estest()*(-t), so estest
%   is called without any arguments. Use "estest -t" or "estest -t *"
%   instead.
%     -d   display results
%     -r   recalculate and store regression data
%     -t   report timings
%     -c   report code coverage
%     -l   list lines of code not covered by the tests (with -c)
%
%   Either the command syntax, as above, or the function syntax, e.g.
%   estest('pepper*','-t'), can be used. A cell array of test names can be
%   given as well, e.g. estest({'pepper_c2h','sop_spinonehalf'}).
%
%   estest runs the tests in the tests folder of the EasySpin source
%   repository and can be called from any folder. To be recognized as a
%   test, a file name must contain an underscore (_). Files without _ are
%   ignored.
%
%   A test function has one of these signatures:
%     ok = test()                    direct test
%     ok = test(opt)                 direct test that responds to options
%     [ok,data] = test(opt,refdata)  regression test with reference data
%
%   ok is true/false, or an array of true/false with one element per
%   subtest. The test passes if all elements are true. An empty ok means
%   the test was not run (e.g. no reference data available).
%
%   opt is a structure with the following fields:
%     opt.Display     true/false - whether the test should plot/print (-d)
%     opt.Regenerate  true/false - whether reference data is being
%                      regenerated (-r)
%     opt.Verbosity   1 with -d, 0 otherwise - can be passed on to EasySpin
%                      functions for additional logging
%
%   For regression tests, data is stored in data/<testname>.mat when
%   estest is called with -r, and passed back as refdata in later runs.
%   With -r, refdata is empty.
%
%   To get the results, request an output: out = estest(...). out.outcomes
%   contains one code per test (0 pass, 1 failed, 2 crashed, 3 not tested),
%   and out.Results contains the details.
%
%   See README.md in the tests folder for more information.

function out = estest(varargin)

% Check whether EasySpin is on the MATLAB path
EasySpinPath = fileparts(which('easyspin'));
if isempty(EasySpinPath)
  error('EasySpin is not on the MATLAB path!');
end

% Change to tests folder, and change back when done (also on error or Ctrl+C)
testsDir = fullfile(fileparts(EasySpinPath),'tests');
if ~isfolder(testsDir)
  error('Test folder %s not found. estest requires the EasySpin source repository.',testsDir);
end
oldDir = cd(testsDir);
restoreDir = onCleanup(@()cd(oldDir));  %#ok<NASGU>

fid = 1;  % output to command window

% Separate options (starting with -) from test name patterns
flags = '';
patterns = {};
for iArg = 1:nargin
  arg = varargin{iArg};
  if isstring(arg)
    arg = cellstr(arg);
  end
  if iscell(arg)
    patterns = [patterns arg(:).'];  %#ok<AGROW>
  elseif ischar(arg)
    if startsWith(arg,'-')
      flags = [flags arg(2:end)];  %#ok<AGROW>
    else
      patterns{end+1} = arg;  %#ok<AGROW>
    end
  else
    error('Inputs must be strings/character arrays, or cell arrays of such.');
  end
end
unknownFlags = setdiff(flags,'drtcl');
if ~isempty(unknownFlags)
  error('Unknown option(s): %s. Valid options are -d, -r, -t, -c, -l.',...
    strjoin(cellstr(unknownFlags(:)),', '));
end
runAll = isempty(patterns) || all(strcmp(patterns,'*'));
if isempty(patterns)
  patterns = {'*'};
end

% Options to pass along to test functions
Opt.Display = any(flags=='d');
Opt.Regenerate = any(flags=='r');
Opt.Verbosity = double(Opt.Display);

displayTimings = any(flags=='t');
runCodeCoverageAnalysis = any(flags=='c');
listMissedLines = any(flags=='l');


% Assemble list of tests to be run: all m files in folder with _ in their
% name that match any of the patterns (exact match unless * is used)
fileList = dir('*.m');
allTestNames = {fileList.name};
allTestNames = allTestNames(contains(allTestNames,'_'));
allTestNames = erase(allTestNames,regexpPattern('\.m$'));
selectedTests = {};
for p = 1:numel(patterns)
  pattern = regexprep(patterns{p},'\.m$','');
  regex = ['^' regexptranslate('wildcard',pattern) '$'];
  matches = allTestNames(~cellfun(@isempty,regexp(allTestNames,regex,'once')));
  if isempty(matches)
    fprintf('No tests matching ''%s''.\n',pattern);
  end
  selectedTests = [selectedTests matches];  %#ok<AGROW>
end
if isempty(selectedTests)
  if nargout==1
    out.Results = struct([]);
    out.outcomes = [];
  end
  return
end
testFileNames = strcat(unique(selectedTests).','.m');

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

% List all EasySpin functions (not including private ones)
if runCodeCoverageAnalysis
  Files = dir(fullfile(EasySpinPath,'*.m'));
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
  startTime = tic;
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
  timeElapsed(iTest) = toc(startTime);

  % Wait for keypress (after timing, so waiting is not included)
  if Opt.Display
    if iTest<numel(testFileNames)
      pause;
    end
  end  
  
  % Retrieve profiler summary and turn profiler off
  if runCodeCoverageAnalysis
    p = profile('info');
    profile off
    
    % Make list of files of all profiled function calls
    executedFiles = {p.FunctionTable(:).FileName};
    % Analyze code coverage of each API function
    for n = 1:length(Files)
      fcnFile = fullfile(Files(n).folder,Files(n).name);
      if ispc
        pos = find(strcmpi(executedFiles,fcnFile));
      else
        pos = find(strcmp(executedFiles,fcnFile));
      end
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

  nBlanks = max(47-length(nameStr),0);
  testFile = strrep(fullfile(testsDir,[nameStr '.m']),'\','/');
  nameStrLink = sprintf('<a href="matlab: edit(''%s'')">%s</a>%s',testFile,nameStr,repmat(' ',1,nBlanks));
  strLink = sprintf('%s  %-12s%-8s%s\n%s',...
       nameStrLink,typeStr,outcomeStr,timeStr,errorStr);

  testResults(iTest).msg = str;
  testResults(iTest).msgLink = strLink;

  fprintf(fid,'%s',strLink);
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
    Code = fileread(fullfile(Path,fcnName));
    % Account for lines missed by callstats
    if covered > runnable
        totalRunnable = totalRunnable + covered;
    else
        totalRunnable = totalRunnable + runnable;
    end
    missed = [];
    if listMissedLines
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
    if runAll || coverage~=0
        fprintf('%-20s%-18s%5.1f%% %18s  %s\n',fcnName,' ',coverage,'Lines missing:',mat2str(missed))
    end
  end
  totalCoverage = totalCovered/totalRunnable*100;
  if runAll
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
failedOrCrashed = find(allOutcomes==1 | allOutcomes==2);
if ~isempty(failedOrCrashed)
  fprintf(fid,'-----------------------------------------------------------------------\n');
  for iTest = failedOrCrashed
    fprintf(fid,'%s',testResults(iTest).msgLink);
  end
end

nPasses = sum(allOutcomes==0);
nFailures = sum(allOutcomes==1);
nCrashes = sum(allOutcomes==2);
fprintf(fid,'-----------------------------------------------------------------------\n');
fprintf(fid,'%d passes, %d failures, %d crashes\n',nPasses,nFailures,nCrashes);
fprintf(fid,'-----------------------------------------------------------------------\n');

% Assemble list of all failed and crashed tests and provide link to rerun them.
if nFailures+nCrashes>0
  testList = cell(1,nFailures+nCrashes);
  idx = 1;
  for iTest = failedOrCrashed
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
