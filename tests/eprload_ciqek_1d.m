function ok = test(opt)   %#ok

% Read CIQTEK spectrometer files with 1D data
%-------------------------------------------------------------------------------

baseDir = 'eprfiles/ciqtek/';
fileNames{1} = 'CIQTEK_1D.epr';
nFiles = numel(fileNames);


readerr = false(1,nFiles);
for iFile = 1:nFiles
  thisFile = fileNames{iFile};
  if opt.Display
    fprintf('file: %s\n',thisFile);
  end
  thisFileName = [baseDir thisFile];

  clear x data pars
  data = eprload(thisFileName);  %#ok
  [x,data] = eprload(thisFileName);   %#ok
  [x,data,pars] = eprload(thisFileName);

  % x must be a 1D array
  if ~isvector(x)
    readerr(iFile) = true;
  end

  % data must be a 1D array
  if ~isvector(data)
    readerr(iFile) = true;
  end
  
  % pars must be a structure
  if ~isstruct(pars)
    readerr(iFile) = true;
  end

  if opt.Display
    fprintf('data size: %d x %d\n',size(x,1),size(x,2));
    plot(x,data);
    pause
  end

end

ok = ~readerr;
