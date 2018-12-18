function [err,data] = test(opt,olddata)

% Read Magnettech files containing 2D sweeps (field + other parameter
%-------------------------------------------------

BaseDir = 'eprfiles/magnettech/';

Files{1} = 'temperatureSweep.xml';
Files{2} = 'powerSweep.xml';
Files{3} = 'modulationSweep.xml';

for iFile = 1:numel(Files)
  fullFileName = [BaseDir Files{iFile}];
  [x,data,pars] = eprload(fullFileName);
  ok(iFile) = true;
end

err = any(~ok);

data = [];
