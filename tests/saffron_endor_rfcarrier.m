function ok = test() %#ok<FNDEF>

% ENDOR measures the echo as a function of the RF carrier frequency
% (Exp.rfFreq, swept in Exp.Dim1). The RF pulse frequency (Pulse.Frequency)
% is an offset from the carrier: an offset f0 shifts the ENDOR spectrum by
% -f0 along the carrier axis. Chirped RF pulses are defined like chirped
% microwave pulses, with a frequency band relative to the carrier.

Sys.Nucs = '1H';
Sys.A = [2 2 5];  % MHz

p90.Flip = pi/2;
rf.Channel = 'rf';
rf.tp = 5;  % µs
rf.Flip = pi;

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz, needed by the full-space engine (chirped RF pulse)
Exp.MolFrame = [0 0.2 0];
Exp.CrystalSymmetry = 1;
Exp.rfFreq = 10;  % MHz
Exp.Sequence = {p90 0.15 p90 0 rf 0 p90 0.15};
Exp.nPoints = 201;
Exp.Dim1 = {'rfFreq',0.05};
Exp.Pathways = '+00-';
Exp.PhaseCycle = {[],[],[0 -1 1; 0 1 0],[]};

[x,y0] = saffron(Sys,Exp);

% carrier offset of the pulse: shifts the spectrum
f0 = 1;  % MHz
E1 = Exp;
E1.Sequence{5}.Frequency = f0;
[~,y1] = saffron(Sys,E1);
shift = round(f0/0.05);
ok(1) = areequal(y1(1:end-shift),y0(1+shift:end),1e-10,'rel');

% chirped RF pulse (adiabatic inversion over 0.4 MHz, full-space engine):
% broader peaks with a plateau, centered at the same frequencies
chirp.Channel = 'rf';
chirp.Type = 'quartersin/linear';
chirp.tp = 20;  % µs
chirp.trise = 4;  % µs
chirp.Frequency = [-0.2 0.2];  % MHz
chirp.Qcrit = 5;
E2 = Exp;
E2.Sequence{5} = chirp;
[~,y2,info] = saffron(Sys,E2);
c0 = sum(x.*abs(y0))/sum(abs(y0));
c2 = sum(x.*abs(y2))/sum(abs(y2));
ok(2) = strcmp(info.Method,'fullspace') && contains(info.MethodReason,'shaped or frequency-swept RF pulses');
ok(3) = abs(c2-c0)<0.1 && sum(abs(y2)>0.5*max(abs(y2))) > sum(abs(y0)>0.5*max(abs(y0)));

% the subspace engine does not support shaped RF pulses
try
  saffron(Sys,E2,struct('Method','subspace'));
  ok(4) = false;
catch err
  ok(4) = contains(err.message,'shaped or frequency-swept RF pulses');
end
