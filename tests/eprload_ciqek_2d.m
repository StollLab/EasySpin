function ok = test(opt)   %#ok

% Read CIQTEK spectrometer files with 2D data
%-------------------------------------------------------------------------------

baseDir = 'eprfiles/ciqtek/';
fileNames{1} = 'CIQTEK_2D.epr';
nFiles = numel(fileNames);

readerr = false(1,nFiles);
for iFile = 1:nFiles
  thisFile = fileNames{iFile};
  if opt.Display
    fprintf('file: %s\n',thisFile);
  end
  thisFileName = [baseDir thisFile];

  clear xy data pars
  data = eprload(thisFileName);  %#ok
  [xy,data] = eprload(thisFileName);   %#ok
  [xy,data,pars] = eprload(thisFileName);

  % x axis must match first axis of 2D data array
  if size(data,1)~=numel(xy{1})
    readerr(iFile) = true;
  end
  % y axis must match second axis of 2D data array
  if size(data,2)~=numel(xy{2})
    readerr(iFile) = true;
  end
  % pars must be a structure
  if ~isstruct(pars)
    readerr(iFile) = true;
  end

  if opt.Display
    fprintf('data size: %d x %d\n',size(xy,1),size(xy,2));
    plot(xy,data);
    pause
  end

end

ok = ~readerr;
