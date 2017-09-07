% Custom function for simulating slow-motion X- and Q-band EPR spectra
function y = chili_multifreq(Sys, Exp, SimOpt)

% X-band spectrum
Sys.lw = Sys.lwX;

Exp.CenterSweep = Exp.CenterSweepX;
Exp.nPoints = Exp.nPointsX;
Exp.mwFreq = Exp.mwFreqX;

y = chili(Sys, Exp, SimOpt);
yX = y/max(y);  % normalize so that spectra can be viewed on same scale 
                % while fitting

% Q-band spectrum
Sys.lw = Sys.lwQ;

Exp.CenterSweep = Exp.CenterSweepQ;
Exp.nPoints = Exp.nPointsQ;
Exp.mwFreq = Exp.mwFreqQ;

y = chili(Sys, Exp, SimOpt);
yQ = y/max(y);

% Combine X- and Q-band spectra for fitting
yXref = Sys.dataX;
yQref = Sys.dataQ;

% rescale each spectrum individually here since, internally, esfit will
% rescale the concatenated output y
y = [rescale(yX,yXref,'lsq'), rescale(yQ,yQref,'lsq')];

return