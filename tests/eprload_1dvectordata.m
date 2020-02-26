function ok = test(opt)

BaseDir = 'eprfiles/';
Files{1} = [BaseDir 'calib_10dB_RT_9370-9440mT_lin_1G.DTA'];
nDataVals(1) = 2;

%-------------------------------------------------
% Read datasets with vector data
% (multiple reals or complex scalars per parameter point
%-------------------------------------------------

for iFile = 1:numel(Files)
  if (opt.Display)
    disp(Files{iFile});
  end
  readerr(iFile) = false;
  data = eprload(Files{iFile});
  [x,data] = eprload(Files{iFile});
  [x,data,pars] = eprload(Files{iFile});
  if ~iscell(data)
    readerr(iFile) = true;
  end
  if numel(data)~=nDataVals(iFile)
    readerr(iFile) = true;
  end
  if (opt.Display)
    for k = 1:nDataVals
      subplot(nDataVals,1,k);
      cla
      if isreal(data{k})
        plot(real(data{k}));
      else
        %x = 1:numel(data{k});
        plot(x,imag(data{k}),'r',x,real(data{k}),'b');
      end
    end
    pause
  end
end

ok = ~any(readerr);
