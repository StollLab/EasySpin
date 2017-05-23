function [err,data] = test(opt,olddata)

% Read Adani spectrometer files
%-------------------------------------------------------------------------------
% JSON format from e-Spinoza software

% Skip test if Matlab version is older than R2016b
% Before that, the Matlab-internal JSON parser does not work.
if verLessThan('matlab','9.1')
  err = false;
  data = [];
  return
end

BaseDir = 'eprfiles/adani/';
Files{1} = '1d.json';
Files{2} = '2d.json';

for iFile = 1:numel(Files)
  thisFile = Files{iFile};
  if (opt.Display)
    disp(thisFile);
  end
  clear data pars;
  thisFileName = [BaseDir thisFile];
 
  data = eprload(thisFileName);
  [x,data] = eprload(thisFileName);
  [x,data,pars] = eprload(thisFileName);
  
  readerr(iFile) = false;
  if ~iscell(data)
    if any(isnan(data(:)))
      readerr(iFile) = true;
      disp([  '   ' Files(iFile).name]);
    end
  end
  
  if (opt.Display)
    plot(x,data);
    pause
  end
  
end

err = readerr;

data = [];
