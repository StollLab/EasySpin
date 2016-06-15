% Impact of different scaling methods on fitted parameters

clear

% If a noisy experimental spectrum is to be fitted, the way
% the simulated spectrum is scaled before comparing it to the
% experimental one is important.

% Let's generate a simple spectrum with some noise.
Sys.g = 2;
Sys.lw = 2;
Exp.mwFreq = 9.5;
Exp.CenterSweep = [339.4 20];
[x,y] = pepper(Sys,Exp);
y = addnoise(y,10,'n');

% If we fit a linewidth to this spectrum, we get different
% results depending on whether we use 'maxabs', 'minmax' or
% 'lsq' fitting.
Sys0 = Sys;
Sys0.lw = 1.6;   % slightly offset the line width from its correct value
Vary.lw = 0.8;

FitOpt.Method = 'simplex int';
FitOpt.PrintLevel = 0;

FitOpt.Scaling = 'minmax';
f1 = esfit('pepper',y,Sys0,Vary,Exp,[],FitOpt);
FitOpt.Scaling = 'maxabs';
f2 = esfit('pepper',y,Sys0,Vary,Exp,[],FitOpt);
FitOpt.Scaling = 'lsq1';
f3 = esfit('pepper',y,Sys0,Vary,Exp,[],FitOpt);

% If you run this code several times, you will find that
% 'lsq' seems to be the best method.
[Sys.lw f1.lw f2.lw f3.lw]
