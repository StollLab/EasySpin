function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
HS.Type = 'sech/uniformQ';
HS.beta = 10;
HS.n = [10 10];

Exp.t = [0.2 0.5 0.2];
Exp.Pulses = {HS 0 HS};
Exp.Flip = [pi/2 pi];
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.1 0.1];
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 

% Options ---------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;
Opt.SimulationMode = 'FrameShift';

% To test I --------------------------
Exp.nPoints = [3];
Exp.Dim1 = {'p1.n(2)' -3};

[~, signal1] = spidyan(Sys,Exp,Opt);

% To test II --------------------------
Exp.Dim1 = {'p1.n(2)' [0 -3 -6]};

[~, signal2] = spidyan(Sys,Exp,Opt);


if ~areequal(signal1,signal2,1e-4)
  err = 1;
else
  err = 0;
end

data = [];