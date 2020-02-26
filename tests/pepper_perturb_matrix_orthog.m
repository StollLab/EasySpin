function ok = test(opt)

% test ful gmatrix spec with perturbation theory
Sys.g = [2 4 6];
Sys.lw = 5;
%Sys.Nucs = '1H';
%Sys.A = [500 400];
Exp.mwFreq = 9.5;
Exp.Range = [80 380];
Exp.Harmonic = 1;

Opt.Method = 'matrix';
[x0,y0]=pepper(Sys,Exp,Opt);
Opt.Method = 'perturb2';
[x1,y1]=pepper(Sys,Exp,Opt);


if opt.Display
  plot(x0,y0,x1,y1);
  legend('matrix','perturb2');
end

ok = areequal(y0,y1,5e-3,'abs');
