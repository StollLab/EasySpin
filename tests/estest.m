% estest    Testing engine for EasySpin
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
%   Run all tests including timings:   estest('*','t')
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
  error('No test functions matching the pattern %s',FileMask);
end

TestFileNames = sort({FileList.name});

fprintf(fid,'=======================================================================\n');
fprintf(fid,'EasySpin test set                      %s\n(MATLAB %s)\n',datestr(now),version);
fprintf(fid,'EasySpin location: %s\n',EasySpinPath);
fprintf(fid,'=======================================================================\n');
fprintf(fid,'Display: %d, Regenerate: %d, Verbosity: %d\n',...
  Opt.Display,Opt.Regenerate, Opt.Verbosity);
fprintf(fid,'-----------------------------------------------------------------------\n');

% Codes for test outcomes:
%    0   test passed
%   +1   test failed
%   +2   test crashed
%   +3   not tested

OutcomeStrings = {'pass','failed','crashed','not tested'};

%get path to Easyspin functions folder
path = fileparts(which('estest'));
path = path(1:end-length('\tests'));

%list all the functions, including private
Files = dir(fullfile(path,'easyspin','*.m'));
ExecutedLines = repmat({[]},length(Files),1);

for iTest = 1:numel(TestFileNames)
  
  thisTestName = TestFileNames{iTest}(1:end-2);

  if Opt.Display
    clf
    set(gcf,'Name',thisTestName);
    drawnow
  end
  
  % Load, or regenerate, comparison data
  olddata = [];
  TestDataFile = ['data/' thisTestName '.mat'];
  if exist(TestDataFile,'file')
    if Opt.Regenerate
      delete(TestDataFile);
      olddata = [];
    else
      try
        olddata = load(TestDataFile,'data');
        olddata = olddata.data;
      catch
        error('Could not load data for test ''%s''.',thisTestName);
      end
    end
  end
  
  % Run test, catch any errors
  testFcn = str2func(thisTestName);
    
  %Clear and start profiling code run by tests
  profile clear
  profile on
  tic
  try
    [err,data] = testFcn(Opt,olddata);
    if Opt.Display
      if iTest<numel(TestFileNames), pause; end
    end
    % if test returns empty err, then treat it as not tested
    if isempty(err)
      err = 3; % not tested
    else
      err = any(err~=0);
    end
    errorInfo = [];
    errorStr = '';
  catch exception
    data = [];
    err = 2;
    errorInfo = exception;
    errorStr = getReport(errorInfo);
    errorStr = ['    ' regexprep(errorStr,'\n','\n    ') char(10)];
  end
  time_used(iTest) = toc;
  
   %retrieve profiler summary
  p = profile('info');
  %and turn it off
  profile off
  
  %Make list of all profiled function calls
  ExecutedFcns = [{p.FunctionTable(:).CompleteName}];
  %Analyze code coverage of each API function
  for n = 1:length(Files)
      FcnName = Files(n).name;
      pos = find(contains(ExecutedFcns,FcnName));
      if ~isempty(pos)
          %initialize containers
          for i=1:length(pos)
              %get executed lines in profiler
              tmp = p.FunctionTable(pos(i)).ExecutedLines;
              container = ExecutedLines{n};
              container(end+1:end+length(tmp(:, 1))) = tmp(:, 1);
              ExecutedLines{n} = container;
          end
      end
  end
  
  isRegressionTest = ~isempty(data);
  saveTestData = isRegressionTest && isempty(olddata);  
  if saveTestData
    save(TestDataFile,'data');
  end
  
  testResults(iTest).err = double(err);
  testResults(iTest).name = thisTestName;
  testResults(iTest).errorData = errorInfo;
  
  outcomeStr = OutcomeStrings{testResults(iTest).err+1};
  
  if ~isempty(data)
    typeStr = 'regression';
  else
    typeStr = 'direct';
  end
  
  if displayTimings
    timeStr = sprintf('%0.3f seconds',time_used(iTest));
  else
    timeStr = [];
  end

  str = sprintf('%-45s  %-12s%-8s%s\n%s',...
       testResults(iTest).name,typeStr,outcomeStr,timeStr,errorStr);
  str(str=='\') = '/';
  
  testResults(iTest).msg = str;
  
  fprintf(fid,str);
end



fprintf(fid,'-----------------------------------------------------------------------\n');
fprintf(fid,'Code Coverage Analysis \n');
fprintf(fid,'-----------------------------------------------------------------------\n');

TotalCovered = 0;
TotalRunnable = 0;
%Analyze code coverage of each API function
for n = 1:length(Files)
    FcnName = Files(n).name;
    Path = Files(n).folder;
    RunnableLines = callstats('file_lines',fullfile(Path,FcnName));
    if isempty(RunnableLines)
        RunnableLines = 0;
    end
    TotalRunnable = TotalRunnable + length(unique(RunnableLines));
    Executed = unique(ExecutedLines{n});
    Covered = length(Executed);
    TotalCovered = TotalCovered + Covered;
    Runnable = length(unique(RunnableLines));
    Code = fileread(FcnName);
    if params =='l'
        Missed = RunnableLines;
       for k=1:length(Executed)
           Missed(RunnableLines==Executed(k)) = NaN;
       end
       Missed(isnan(Missed)) = [];
    end
    %Account for unreachable end-statements or lines
    MissedEnds = length(strfind(Code,'error')) + length(strfind(Code,'return')) ...
        + length(strfind(Code,'break')) + length(strfind(Code,'fprintf')) + length(strfind(Code,'plot'));
    if Runnable - Covered <= MissedEnds
        Covered = Runnable;
        Missed = [];
    end
    Coverage = 100*Covered/Runnable;
    %Print to console
    if (~isempty(TestName) && Coverage~=0) || isempty(TestName)
        if params =='l'
            fprintf('%-20s%-18s%3.2f%% %18s  %s \n',FcnName,' ',Coverage,'Lines missing:',mat2str(Missed))
        else
            fprintf('%-20s%-18s%3.2f%%\n',FcnName,' ',Coverage)
        end
    end
end
TotalCoverage = TotalCovered/TotalRunnable*100;

allErrors = [testResults.err];

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
if any(allErrors==1) || any(allErrors==2)
  fprintf(fid,'-----------------------------------------------------------------------\n');
  for iTest = find(allErrors)
    fprintf(fid,testResults(iTest).msg);
  end
end

fprintf(fid,'-----------------------------------------------------------------------\n');
msg = sprintf('%d passes, %d failures, %d crashes\n',sum(allErrors==0),sum(allErrors==1),sum(allErrors==2));
fprintf(fid,msg);
if isempty(TestName)
    fprintf('Total code coverage: %3.2f%%\n',TotalCoverage)
end
fprintf(fid,'-----------------------------------------------------------------------\n');

% Return output if desired
if nargout==1
  out.Results = testResults;
  out.outcomes = allErrors;
end

return
