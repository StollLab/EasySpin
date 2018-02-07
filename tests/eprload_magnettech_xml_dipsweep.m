function [err,y] = test(opt,olddata)

% Read capture of dip sweep (raw data) in new XML format
%--------------------------------------------------------

thisFileName = 'eprfiles/magnettech/DipSweep.xml';
[nu,y,pars] = eprload(thisFileName);
  
err = ~isnumeric(y) || any(isnan(y));

dataloaded = ...
  isfield(pars,'Curves') && ...
  isfield(pars.Curves,'Frequency') && ...
  isfield(pars.Curves,'ADC_24bit');

err = err || ~dataloaded;

if opt.Display
  plot(nu,y);
  axis tight;
  xlabel('frequency (GHz)');
end

y = [];
