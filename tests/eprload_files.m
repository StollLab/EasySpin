function ok = test(opt)

BaseDir = 'eprfiles/';

%-------------------------------------------------
% Test: Read several formats
%-------------------------------------------------

Extensions = {...
'*.dsc','*.dta','*.DSC','*.DTA',...
'*.spc','*.par','*.PAR','*.SPC','*.eco',...
'*.d00','*.d01'...
};

Files = [];
for k = 1:numel(Extensions)
  Files = [Files;dir([BaseDir Extensions{k}])];
end

for iFile = 1:numel(Files)
  if (opt.Display)
    disp(Files(iFile).name);
  end
  clear data pars;
  data = eprload([BaseDir Files(iFile).name]);
  [x,data] = eprload([BaseDir Files(iFile).name]);
  [x,data,pars] = eprload([BaseDir Files(iFile).name]);
  readerr(iFile) = false;
  if ~iscell(data)
    if any(isnan(data(:)))
      readerr(iFile) = true;
      disp([  '   ' Files(iFile).name]);
    end
  end
  if (opt.Display)
    subplot(2,1,1);
    plot(real(data));
    subplot(2,1,2);
    plot(imag(data+0i));
    pause
  end
end

ok = ~readerr;
