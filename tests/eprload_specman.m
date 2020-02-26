function ok = test(opt)

%-------------------------------------------------
% Read specman files
%-------------------------------------------------

BaseDir = 'eprfiles/specman/';
Files = dir(fullfile(BaseDir,'*.d01'));

for iFile = 1:numel(Files)
  fname = Files(iFile).name;
  fullname = [BaseDir fname];
  if strcmp(fname,'.'), continue; end
  if strcmp(fname,'..'), continue; end
  if (opt.Display), disp(fname); end
  
  data = eprload(fullname);
  [x,data] = eprload(fullname);
  [x,data,pars] = eprload(fullname);
  
  readerr(iFile) = false;
  if ~iscell(data)
    if any(isnan(data(:)))
      readerr(iFile) = true;
      disp([  '   ' Files(iFile).name]);
    end
  end
  
  if opt.Display
    subplot(2,1,1);
    plot(x,real(data));
    subplot(2,1,2);
    plot(x,imag(data+0i));
    pause
  end
end

ok = ~readerr;
