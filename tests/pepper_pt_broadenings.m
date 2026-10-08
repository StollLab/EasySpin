function ok = test(opt)

% Test whether H strain works equally
% with matrix diagonalization and perturbation theory
%-------------------------------------------------------------

Exp.mwFreq = 9.5;
Exp.Range = [200 400];
Sys.Nucs = '1H';
Sys.g = [2 2.3 2.7];
Sys.A = [350 400];
Sys.HStrain = [200 250 300];
Opt.Method = 'matrix';
[x,y0] = pepper(Sys,Exp,Opt);
Opt.Method = 'perturb2';
[x,y2] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(x,y0,x,y2,'g')
  legend('matrix','perturb2');
  legend boxoff
end

ok = areequal(y0,y2,5e-3,'abs');
