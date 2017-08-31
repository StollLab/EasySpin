% estest    Testing engine for EasySpin
%
%   Usage:
%     estest            run all tests
%     estest adsf       run all tests whose name starts with asdf
%     estest asdf d     run all tests whose name starts with asdf and
%                       display results
%     estest asdf r     evaluate all tests whose name starts with asdf and
%                       recalculate and store regression data
%     estest asdf t     evaluate all tests whose name starts with asdf and
%                       report timings
%
%   Either the command syntax as above or the function syntax, e.g.
%   estest('asdf','t'), can be used.
%
%   Run all tests including timings:   estest('*','t')
%
%   All test files must have an underscore _ in their filename.

function out = estest(TestName,params)

% Check whether EasySpin is on the Matlab path
EasySpinPath = fileparts(which('easyspin'));
if isempty(EasySpinPath)
  error('EasySpin is not on the Matlab path!');
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
fprintf(fid,'EasySpin test set                      %s\n(Matlab %s)\n',datestr(now),version);
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

for iTest = 1:numel(TestFileNames)
  
  if (Opt.Display)
    clf; drawnow;
  end

  thisTest = TestFileNames{iTest}(1:end-2);
  
  olddata = [];
  TestDataFile = ['data/' thisTest '.mat'];
  if exist(TestDataFile,'file')
    if Opt.Regenerate
      delete(TestDataFile);
    else
      try
        olddata = load(TestDataFile,'data');
        olddata = olddata.data;
      catch
        lasterr
        olddata = [];
      end
    end
  end

  lasterr('');
  tic
  try
    [err,data] = feval(thisTest,Opt,olddata);
    if (Opt.Display)
      if (iTest<numel(TestFileNames)), pause; end
    end
    % if test returns empty err, then treat it as not tested
    if isempty(err)
      err = 3; % not tested
    else
      err = any(err~=0);
    end
  catch
    data = [];
    err = 2;
  end
  time_used(iTest) = toc;
  testError = lasterr;
  
  saveTestData = ~isempty(data) && (isempty(olddata) || Opt.Regenerate);
  
  if saveTestData
    if verLessThan('matlab','7')
      save(TestDataFile,'data');
    else
      save(TestDataFile,'data','-V6');
    end
  end
    
  Results(iTest).err = double(err);
  Results(iTest).name = thisTest;
  Results(iTest).errmsg = testError;
  if ~isempty(testError)
    errStr = sprintf('%s\n',testError);
  else
    errStr = '';
  end
  resultStr = OutcomeStrings{Results(iTest).err+1};
  
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
  
  str = sprintf('%-36s  %-12s%-8s%s\n%s',...
       Results(iTest).name,typeStr,resultStr,timeStr,errStr);
  str(str=='\') = '/';
  
  Results(iTest).msg = str;
  
  fprintf(fid,str);
end

allErrors = [Results.err];

% Display timings of slowest tests
if displayTimings
  fprintf(fid,'-----------------------------------------------------------------------\n');
  fprintf(fid,'Total test time:                        %7.3f seconds\n',sum(time_used));
  fprintf(fid,'Slowest tests:\n');
  [time,iTest] = sort(time_used,'descend');
  for q = 1:min(10,numel(time))
    fprintf(fid,'%-36s    %7.3f seconds\n',Results(iTest(q)).name,time(q));
  end
end

% Display all tests that failed
if any(allErrors==1) || any(allErrors==2)
  fprintf(fid,'-----------------------------------------------------------------------\n');
  for iTest = find(allErrors)
    fprintf(fid,Results(iTest).msg);
  end
end

fprintf(fid,'-----------------------------------------------------------------------\n');
msg = sprintf('%d passes, %d failures, %d crashes\n',sum(allErrors==0),sum(allErrors==1),sum(allErrors==2));
fprintf(fid,msg);
fprintf(fid,'-----------------------------------------------------------------------\n');

% Return output if desired
if nargout==1
  out.Results = Results;
  out.outcomes = allErrors;
end

return
