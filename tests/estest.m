% estest    Testing engine for EasySpin

function estest(TestName,params)

% Check whether EasySpin is on the Matlab path
Espath = fileparts(which('sop'));
if isempty(Espath)
  error('EasySpin not on the Matlab path!');
end

if (nargin<2)
  params = '.';
end

Opt.Display = any(params=='d');
Opt.Regenerate = any(params=='r');
Opt.Verbosity = Opt.Display;

OldPath = path;
OldDir = pwd;

LogLevel = 1;

if (LogLevel>0)
  fprintf('Display: %d, Regenerate: %d, Verbosity: %d\n',...
    Opt.Display,Opt.Regenerate, Opt.Verbosity);
end

fid = 1; % standard output

if (nargin==0), TestName = '*'; end

if strfind(TestName,'_')
  FileMask = [TestName '*.m'];
else
  FileMask = [TestName '*_*.m'];
end

FileList = dir(FileMask);

%for f=numel(FileList):-1:1
%  if ~isempty(strfind(FileList(f).name,'_test')), FileList(f) = []; end
%end

if numel(FileList)==0
  error('No test functions matching the pattern %s',FileMask);
end

for f = 1:numel(FileList)
  TestNames{f} = FileList(f).name;
end

[unused,idx] = sort({FileList.name});
TestNames = TestNames(idx);


if (LogLevel>0)
  fprintf(fid,'=======================================================================\n');
  fprintf(fid,'EasySpin test suite                    %s\n(Matlab %s)\n',datestr(now),version);
  fprintf(fid,'EasySpin location: %s\n',Espath);
  fprintf(fid,'=======================================================================\n');
end

% Error codes:
%  0   test passed
% +1   test failed
% +2   test crashed
% +3   not tested

ErrorStrings = {'ok','failed','crashed','not tested'};

for iTest = 1:numel(TestNames)
  
  if (Opt.Display)
    clf; drawnow;
  end

  thisTest = TestNames{iTest}(1:end-2);
  
  olddata = [];
  TestDataFile = ['data/' thisTest '.mat'];
  if exist(TestDataFile,'file')
    if (Opt.Regenerate)
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
      if (iTest<numel(TestNames)), pause; end
    end
    if isempty(err)
      err = 3;
    else
      err = any(err~=0);
    end
  catch
    data = [];
    err = 2;
  end
  time_used(iTest) = toc;
  
  SaveTestData = ~isempty(data) & (isempty(olddata) | Opt.Regenerate);
  
  if (SaveTestData)
    vv = version;
    if (vv(1)>='7')
      save(TestDataFile,'data','-V6');
    else
      save(TestDataFile,'data');
    end
  end
    
  Results(iTest).err = double(err);
  Results(iTest).name = thisTest;
  Results(iTest).errmsg = lasterr;
  if ~isempty(lasterr)
    errStr = [lasterr sprintf('\n')];
  else
    errStr = '';
  end
  resultStr = ErrorStrings{Results(iTest).err+1};
  if ~isempty(data), typeStr = 'regression'; else typeStr = 'direct'; end
  Results(iTest).msg = [sprintf('%-36s  %-14s%s\n',...
     Results(iTest).name,typeStr,resultStr) errStr];
  
  Results(iTest).msg(Results(iTest).msg=='\') = '/';
  fprintf(fid,Results(iTest).msg);
end
fprintf(fid,'-----------------------------------------------------------------------\n');

allErrors = [Results.err];

msg = sprintf('%d successful tests, %d failures, %d crashes\n',sum(allErrors==0),sum(allErrors==1),sum(allErrors==2));
fprintf(fid,msg);

if sum(allErrors)>0
  fprintf(fid,'-----------------------------------------------------------------------\n');
  for iTest = find(allErrors)
    fprintf(fid,Results(iTest).msg);
  end
end

fprintf(fid,'-----------------------------------------------------------------------\n');

path(OldPath);
cd(OldDir);

return
