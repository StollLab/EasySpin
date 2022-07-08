% Genetic algorithm combined with local simplex search
%==========================================================================

clear, close

% This examples shows how you can perform a two-stage
% global/local fitting. This "hybrid" method works better
% than either global or local search alone.

% Let's generate a noisy spectrum, acquired at 95 GHz
Sys.g = [2.01 2.04 2.07];
Sys.Nucs = '1H';
Sys.A = [150 90 130];
Sys.lw = 2;
Exp.mwFreq = 95; Exp.CenterSweep = [3330 240];
[x,y] = pepper(Sys,Exp);
y = addnoise(y,200,'n');

% For the starting parameter set of the fitting, we
% modify the g values used for the simulation by small random
% amounts and let the fitting procedure vary them considerably
Sys0 = Sys;
Sys0.g = Sys.g + rand(1,3)*0.005;
Vary.g = [1 1 1]*0.02;

% We perform a global search using the genetic algorithm on
% the integral of the spectrum.
% This can then be followed by a local search using the simplex 
% algorithm, starting from the current best fit 
% (plotted in green).
SimOpt.Method = 'perturb';
FitOpt.Method = 'genetic int';
esfit(y,@pepper,{Sys0,Exp,SimOpt},{Vary},FitOpt);
