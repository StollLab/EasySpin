function ok = test() %#ok<FNDEF>

% An indirect dimension over one element of a pulse parameter (sech pulse
% truncation n(2)): a linear increment and a list of values give the same
% result (single orientation, Sz detected over the whole sequence)

% System ------------------------
Sys.S = 1/2;

% Experiment -------------------
HS.Type = 'sech/uniformQ';
HS.beta = 10;
HS.n = [10 10];
HS.tp = 0.2;
HS.Frequency = 1000* [-0.1 0.1];

HS1 = HS;
HS1.Flip = pi/2;

HS2 = HS;
HS2.Flip = pi;


Exp.Sequence = {HS1 0.5 HS2};
Exp.mwFreq = 33.5;
Exp.Field = mhz2mt(33500);  % mT, resonance field for g = gfree
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetWindow = 'sequence';

% Options ---------------------------
Exp.DetOperator = {'z1'};
% simulation frame at the microwave frequency: long time steps
Opt.SimFreq = 33.5;

% To test I --------------------------
Exp.nPoints = 3;
Exp.Dim1 = {'p1.n(2)' -3};

[~, signal1] = saffron(Sys,Exp,Opt);

% To test II --------------------------
Exp.Dim1 = {'p1.n(2)' [0 -3 -6]};

[~, signal2] = saffron(Sys,Exp,Opt);

ok = areequal(signal1,signal2,1e-4,'abs');
