function ok = test(opt)

% Checks whether product-rule simulations are independent
% of the order of the nuclei in Sys.Nucs.

Exp.Sequence = '3pESEEM';
Exp.Field = 340.2; % mT
Exp.ExciteWidth = 30; % MHz
Exp.mwFreq = 9.738; % GHz

Exp.dt = 0.016; % us
Exp.tau = 0.224; % us
Exp.T = 0.080; % us
Exp.nPoints = 500;
Exp.SampleFrame = [1.1234 2.534 0.5];

Opt.GridSize = 101;
Opt.TimeDomain = 1;
Opt.ProductRule = 1;

A_63Cu = [96 96 570];
A_2H = [-0.553 -0.553 1.128];
Sys1.g = [2.046 2.046 2.279];
Sys1.Nucs = '2H,63Cu';  %WORKS FINE
Sys1.A = [A_2H; A_63Cu]; % MHz

[x1, y1] = saffron(Sys1,Exp,Opt);

Sys2.g = [2.046 2.046 2.279];
Sys2.Nucs = '63Cu,2H'; %DOES NOT WORK
Sys2.A = [A_63Cu; A_2H]; % MHz

[x2, y2] = saffron(Sys2,Exp,Opt);

if opt.Display
  subplot(4,1,[1 2 3]);
  plot(x1,real(y1),x2,real(y2));
  legend(Sys1.Nucs,Sys2.Nucs);
  legend boxoff
  axis tight
  ylabel('echo amplitude');
  subplot(4,1,4);
  plot(x1,real(y1)-real(y2));
  xlabel('time (\mus)');
  axis tight
end

ok = areequal(y1,y2,max(y1)*1e-6,'abs');
