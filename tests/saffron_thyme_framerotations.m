function [err,data] = test(opt,olddata)

angles = [pi/2 pi/4 0];

g = gfree+[0 0.005 0.01];
A = [5 2 1];

Sys.S = [1/2 1/2];
Sys.J = 0;
Sys.Nucs = '1H,1H';

Sys.g = [g; g];
Sys.A = [A 0 0 0; 0 0 0 A];

Sys1 = Sys;
Sys1.gFrame = [angles; angles];
Sys1.AFrame = [angles angles; angles angles];

Sys2 = Sys;

pulse.tp = 0.01;
pulse.Flip = pi;

Exp.Sequence = {pulse};

Exp.DetWindow = [0 0.3];
Exp.mwFreq = 9.3;
Exp.Field = 332;
Exp.DetPhase = 0;

Exp.CrystalOrientation = [0 0 0];

[x1, y1] = saffron(Sys1,Exp);

Exp.CrystalOrientation = angles;

[x2, y2] = saffron(Sys2,Exp);

if (opt.Display)
    clf
    plot(x1,real(y1),'r',x2,real(y2),'b');
    legend({'external rotation' 'roation with saffron'})
    axis tight
    ylabel('Signal')
    xlabel('time [us]');
end

err = ~isequal(y1, y2);

data = [];