function [err,data] = test(opt,olddata)

% Read Bruker DTA/DSC files containing AWG data
% (DSC files are unusually large)
%-------------------------------------------------
BaseDir = 'eprfiles/';

Files{1} = '010_cutpp_10kfs.DTA';


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
