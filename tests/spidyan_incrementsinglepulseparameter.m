function [err,data] = test(opt,olddata)

% System ------------------------
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;

% Experiment -------------------
HS.Type = 'sech/uniformQ';
HS.beta = 10;
HS.n = [10 10];
HS.tp = 0.2;
HS.Frequency = [-0.1 0.1];

HS1 = HS;
HS1.Flip = pi/2;

HS2 = HS;
HS2.Flip = pi;


Exp.Sequence = {HS1 0.5 HS2};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5;
Exp.DetSequence = 1; 

% Options ---------------------------
Exp.DetOperator = {'z1'};
Opt.SimFrequency = 32;

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