function ok = test(opt)

% Read Magnettech spectrometer files (XML format)

BaseDir = 'eprfiles/magnettech/';

Files{1} = 'Mangan_Ausgangslage.xml';

for iFile = 1:numel(Files)
  thisFile = Files{iFile};
  if (opt.Display)
    disp(thisFile);
  end
  clear data pars
  thisFileName = [BaseDir thisFile];
 
  data = eprload(thisFileName);
  [B,data] = eprload(thisFileName);
  [B,data,pars] = eprload(thisFileName);

  ok(iFile) = numel(B)==numel(data) && numel(B)>1 && isfield(pars,'MwFreq');
  if ~iscell(data)
    if any(isnan(data(:)))
      ok(iFile) = true;
      disp([  '   ' Files(iFile).name]);
    end
  end
  if opt.Display
    plot(B,data);
    pause
  end
end
