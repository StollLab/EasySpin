function ok = test(opt)

% Make sure saffron takes Exp.SampleFrame into account

angles = [pi/3 pi/5 0];

g = gfree + [0 0.005 0.01];
A = [5 2 1];

Sys.S = [1/2];
Sys.J = 0;
Sys.Nucs = '1H';

Sys.g = g;
Sys.A = A;

% System with tensor tilts
Sys1 = Sys;
Sys1.gFrame = angles;
Sys1.AFrame = angles;

% System without tensor tilts
Sys2 = Sys;

pulse.tp = 0.01;
pulse.Flip = pi;

Exp.Sequence = {pulse};

Exp.DetWindow = [0 0.3];
Exp.mwFreq = 9.3;
Exp.Field = 332;
Exp.DetPhase = 0;

Exp.MolFrame = [0 0 0];

% System with tensor tilts, but without sample tilt
Exp.SampleFrame = [0 0 0];
[x1, y1] = saffron(Sys1,Exp);

% System without tensor tilts, but with sample tilt
Exp.SampleFrame = angles;
[x2, y2] = saffron(Sys2,Exp);

if opt.Display
  clf
  plot(x1,real(y1),x2,real(y2),'--');
  legend({'tensor tilts' 'smaple tilt'})
  axis tight
  ylabel('Signal')
  xlabel('time (µs)');
end

ok = areequal(y1, y2, 1e-12, 'abs');
