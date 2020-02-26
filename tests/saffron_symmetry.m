function ok = test(opt)

% Check whether symmetry is correctly handled

Exp.Sequence = '3pESEEM';
Exp.Field = 324.9;
Exp.dt = 0.008;
Exp.nPoints = 301;
Exp.tau = 0.001;

Sys.Nucs = '1H';
Sys.A = [3 3 12];

Sys.AFrame = [0 0 0];
[x,y1] = saffron(Sys,Exp);
Sys.AFrame = [167 -78 -47]*pi/180;
[x,y2] = saffron(Sys,Exp);

if opt.Display
  xlabel('time [us]');
  plot(x,real(y1),'r',x,real(y2),'b');
  legend('non tilted','tilted');
  title(mfilename)
end

ok = areequal(y1,y2,1e-3,'rel');
