function ok = test()

% Fast mode: pulse phases act on finite pulses (Pulse.tp) in the same way
% as on ideal pulses, so the ratio of the two signals does not depend on
% the pulse phases

Sys.S = 1/2;
Sys.g = gfree;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz

Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
Exp.ExciteWidth = Inf;
Exp.nPoints = 2;
Exp.Dim1 = {'d1,d2' 0.005};

Opt.GridSize = 7;
Opt.SimulationMode = 'fast';

P90.Flip = pi/2;
P180.Flip = pi;
P90real = P90;
P90real.Type = 'rectangular';
P90real.tp = 0.001;  % µs
P180real = P90real;
P180real.Flip = pi;

phases = [0 0; pi/2 0; pi 0; 0 pi/2; pi/3 pi/5];
ratio = zeros(1,size(phases,1));
for k = 1:size(phases,1)
  P90.Phase = phases(k,1); P180.Phase = phases(k,2);
  P90real.Phase = phases(k,1); P180real.Phase = phases(k,2);
  Exp.Sequence = {P90 0.2 P180 0.2};
  yIdeal = saffron(Sys,Exp,Opt);
  Exp.Sequence = {P90real 0.2 P180real 0.2};
  yReal = saffron(Sys,Exp,Opt);
  ratio(k) = yReal(1)/yIdeal(1);
end

ok = areequal(ratio,ratio(1)*ones(size(ratio)),1e-8,'rel');
