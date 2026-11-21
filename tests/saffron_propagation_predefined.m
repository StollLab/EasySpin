function ok = test() %#ok<FNDEF>

% Predefined experiments with the full method: sizes, rejected fields,
% echo amplitudes from the selected echo pathways, 4pESEEM as the diagonal
% of HYSCORE, and real echoes off resonance

Sys.S = 1/2;
Sys.g = 2;

Exp.Field = 350;
Exp.mwFreq = 9.8;
Exp.tau = 0.1;
Exp.dt = 0.02;

Opt.GridSize = 7;
Opt.Method = 'full';

% HYSCORE with a single number of points: used for both dimensions
Exp.Sequence = 'HYSCORE';
Exp.t1 = 0.1;
Exp.t2 = 0.1;
Exp.nPoints = 4;
[x,y] = saffron(Sys,Exp,Opt);
ok(1) = isequal(size(y),[4 4]) && numel(x{1})==4 && numel(x{2})==4;

% Two-pulse ESEEM
Exp = rmfield(Exp,{'t1','t2'});
Exp.Sequence = '2pESEEM';
[~,y] = saffron(Sys,Exp,Opt);
ok(2) = numel(y)==4;

% Exp.DetSequence is no longer supported
Exp.DetSequence = [0 1];
try
  saffron(Sys,Exp,Opt);
  ok(3) = false;
catch err
  ok(3) = contains(err.message,'DetSequence');
end
Exp = rmfield(Exp,'DetSequence');

% Exp.PhaseCycle cannot be used with predefined experiments
Exp.PhaseCycle = [0 0 0; pi 0 pi];
try
  saffron(Sys,Exp,Opt);
  ok(4) = false;
catch err
  ok(4) = contains(err.message,'Exp.PhaseCycle cannot be used with predefined experiments');
end

% Only the echo pathways contribute: without modulation, the echo is
% 1/2 for 2pESEEM and 1/4 for 3pESEEM and HYSCORE (single orientation, on
% resonance)
Sys.g = gfree;
clear Exp
Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
Exp.MolFrame = [0 0 0];
Exp.SampleFrame = [0 0 0];
Exp.tau = 0.1;  % µs
Exp.dt = 0.02;  % µs
Exp.Sequence = '2pESEEM';
Exp.nPoints = 2;
[~,y2] = saffron(Sys,Exp,Opt);
Exp.Sequence = '3pESEEM';
[~,y3] = saffron(Sys,Exp,Opt);
Exp.Sequence = '4pESEEM';
[~,y4] = saffron(Sys,Exp,Opt);
Exp.Sequence = 'HYSCORE';
Exp.nPoints = [2 2];
[~,yH] = saffron(Sys,Exp,Opt);
ok(5) = areequal(y2(:),1/2*ones(2,1),1e-3,'rel');
ok(6) = areequal(y3(:),1/4*ones(2,1),1e-3,'rel');
ok(7) = areequal(yH(:),1/4*ones(4,1),1e-3,'rel');
ok(8) = areequal(y4(:),1/4*ones(2,1),1e-3,'rel');

% 4pESEEM is the diagonal of HYSCORE (t1 = t2 = T)
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz
clear Exp
Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
Exp.tau = 0.1;  % µs
Exp.dt = 0.02;  % µs
Exp.Sequence = '4pESEEM';
Exp.T = 0.05;  % µs
Exp.nPoints = 3;
[x4,y4] = saffron(Sys,Exp,Opt);
Exp = rmfield(Exp,'T');
Exp.Sequence = 'HYSCORE';
Exp.t1 = 0.05;  % µs
Exp.t2 = 0.05;  % µs
Exp.nPoints = [3 3];
[~,yH] = saffron(Sys,Exp,Opt);
ok(9) = areequal(y4(:),diag(yH),1e-10,'rel');
ok(10) = areequal(x4,[0.05 0.07 0.09],1e-12,'abs');

% Ideal pulses refocus the echo exactly at the detection point: off
% resonance, the echo phase is zero
Sys = rmfield(Sys,{'Nucs','A'});
clear Exp
Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq-0.005,'GHz->mT',gfree);  % 5 MHz off resonance
Exp.MolFrame = [0 0 0];
Exp.SampleFrame = [0 0 0];
Exp.tau = 0.1;  % µs
Exp.dt = 0.02;  % µs
Exp.nPoints = 2;
experiments = {'2pESEEM','3pESEEM','4pESEEM','HYSCORE'};
for k = 1:numel(experiments)
  Exp.Sequence = experiments{k};
  [~,y] = saffron(Sys,Exp,Opt);
  ok(10+k) = all(abs(angle(y(:)))<0.01);
end
