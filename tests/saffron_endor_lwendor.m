function ok = test() %#ok<FNDEF>

% Sys.lwEndor broadens ENDOR spectra by convolution along the RF frequency,
% with both its Gaussian and its Lorentzian width, for all ENDOR
% experiments: predefined and user-defined, with RF carrier sweeps and
% sweeps of the RF pulse frequency offset, with several detection
% operators, and in 2D experiments.

Sys.Nucs = '1H';
Sys.A = 6;  % MHz
Opt.GridSize = 7;
lw = [0.3 0.3];  % MHz

% Predefined Mims ENDOR (subspace engine). The broadening keeps the ends of
% the RF range (the echo without RF pulse) and the area of the ENDOR lines
Exp.Field = 350;  % mT
Exp.Sequence = 'MimsENDOR';
Exp.tau = 0.1;  % µs
Exp.rfFreq = [8 22];  % MHz
Exp.nPoints = 351;
y0 = saffron(Sys,Exp,Opt);
Sys.lwEndor = [0 lw(2)];
yL = saffron(Sys,Exp,Opt);
Sys.lwEndor = lw;
yGL = saffron(Sys,Exp,Opt);
range0 = max(real(y0))-min(real(y0));
area = @(y) sum(real(y(1)-y));
ok(1) = max(abs(yL-y0)) > 0.3*range0;  % Lorentzian width alone has an effect
ok(2) = abs(area(yGL)-area(y0)) < 0.05*area(y0);
ok(3) = max(abs(yGL([1 end])-y0([1 end]))) < 1e-3*range0;

% User-defined Mims ENDOR (full-space engine), with a short RF pulse
mwFreq = 350*gfree*bmagn/planck/1e9/1e3;  % GHz
p90.Flip = pi/2;
RF.Channel = 'rf';
RF.Flip = pi;
RF.tp = 2;  % µs
Exp = struct;
Exp.Field = 350;  % mT
Exp.mwFreq = mwFreq;
Exp.rfFreq = 10;  % MHz
Exp.Sequence = {p90 0.1 p90 0 RF 0 p90 0.1};
Exp.nPoints = 45;
Exp.Dim1 = {'rfFreq' 0.25};
Exp.Pathways = '+00-';
Exp.DetPhase = -pi/2;
Opt.Method = 'fullspace';
Sys.lwEndor = 0;
y0 = saffron(Sys,Exp,Opt);
Sys.lwEndor = lw;
y = saffron(Sys,Exp,Opt);
range0 = max(real(y0))-min(real(y0));
ok(4) = max(abs(y-y0)) > 0.3*range0;

% Sweep of the RF pulse frequency offset instead of the RF carrier: the
% same spectrum
ExpOffset = Exp;
RFoffset = RF;
RFoffset.Frequency = 0;  % MHz
RFoffset.Type = 'rectangular';
ExpOffset.Sequence{5} = RFoffset;
ExpOffset.Dim1 = {'p3.Frequency' 0.25};
yOffset = saffron(Sys,ExpOffset,Opt);
ok(5) = max(abs(yOffset-y)) < 0.02*range0;

% Several detection operators: each one is broadened
ExpOps = Exp;
ExpOps.DetOperator = {'+1','z1'};
yOps = saffron(Sys,ExpOps,Opt);
ExpOps.DetOperator = '+1';
yOp1 = saffron(Sys,ExpOps,Opt);
ok(6) = areequal(yOps(:,1).',yOp1,1e-10,'rel');

% 2D, RF carrier x tau: each RF trace is broadened
Exp2D = Exp;
Exp2D.nPoints = [45 2];
Exp2D.Dim2 = {'d1,d4' 0.03};
y2D = saffron(Sys,Exp2D,Opt);
ok(7) = areequal(y2D(:,1).',y,1e-10,'rel');

% RF frequencies as a list of values (subspace engine): the same as with a
% step if equally spaced, an error otherwise
Opt.Method = 'subspace';
Exp = rmfield(Exp,'mwFreq');
yStep = saffron(Sys,Exp,Opt);
Exp.Dim1 = {'rfFreq' (0:44)*0.25};
yList = saffron(Sys,Exp,Opt);
ok(8) = areequal(yList,yStep,1e-10,'rel');
Exp.Dim1{2}(end) = 11.2;
try
  saffron(Sys,Exp,Opt);
  ok(9) = false;
catch me
  ok(9) = contains(me.message,'equally spaced');
end

% all points at the same RF frequency: no broadening
Exp.nPoints = 5;
Exp.Dim1 = {'rfFreq' 0};
yZero = saffron(Sys,Exp,Opt);
Sys.lwEndor = 0;
yZero0 = saffron(Sys,Exp,Opt);
ok(10) = areequal(yZero,yZero0,1e-12,'rel');
