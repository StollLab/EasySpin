function ok = test(opt)

% saffron support for full hyperfine matrices

Exp.Sequence = '3pESEEM';
Exp.tau = 0.001;
Exp.dt = 0.010;
Exp.Field = 350;
Exp.Range =  [0 30];
Opt.GridSymmetry = 'Ci';

Sys.Nucs = '1H';
Sys.lwEndor = 0.1;

Adiag = [2 5 10];
AFrame = [20 60 0]*pi/180;
Ra = erot(AFrame);
A = Ra.'*diag(Adiag)*Ra;

Sys.A = A;
[x,y1]=saffron(Sys,Exp,Opt);

Sys.A = Adiag;
Sys.AFrame = AFrame;
[x,y2]=saffron(Sys,Exp,Opt);

y1=real(y1);
y2=real(y2);

if (opt.Display)
  plot(x,y1,'r',x,y2,'b');
  xlabel('time [us]');
  legend('diag,angle','full matrix');
  title(mfilename)
end

ok = areequal(y1,y2,1e-10,'rel');
