function [err,data] = test(opt,olddata)

% Read Magnettech spectrometer files (new XML format)
%-------------------------------------------------

BaseDir = 'eprfiles/magnettech/';

FileName = 'MINST_050_SK1087.xml';

fullFileName = [BaseDir FileName];
[x,data,pars] = eprload(fullFileName);

ok = isfield(pars,'Curves') && isfield(pars.Curves,'MW_AbsorptionSinus') &&...
  isfield(pars.Curves,'MW_AbsorptionCosinus');
err = ~ok;

if (opt.Display)
  s = pars.Curves.MW_AbsorptionSinus;
  c = pars.Curves.MW_AbsorptionCosinus;
  a = pars.Curves.MW_Absorption;
  plot(s.t,s.data,c.t,c.data,a.t,a.data);
  axis tight
  legend('MW\_AbsorptionSinus','MW\_AbsorptionCosinus','MW\_Absorption');
  legend boxoff
end

data = [];
