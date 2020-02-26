function ok = test(opt)

% full Q tensor, two nuclei

Qpv1 = [-1 -1 2]*0.1;
Qpv2 = [-1.1 -0.9 2]*0.08;
QFrame1 = [10 30 10]*pi/180;
QFrame2 = [30 -40 90]*pi/180;
RQ1 = erot(QFrame1);
RQ2 = erot(QFrame2);

Sys1.Nucs = '14N,14N';
Sys1.A = [1 2.5];
Sys1.Q = [RQ1.'*diag(Qpv1)*RQ1; RQ2.'*diag(Qpv2)*RQ2];

Sys2.Nucs = '14N,14N';
Sys2.A = [1 2.5];
Sys2.Q = [Qpv1; Qpv2];
Sys2.QFrame = [QFrame1; QFrame2];

Sys1.lwEndor = 0.05;
Sys2.lwEndor = 0.05;

Exp.Range = [0 7];
Exp.Field = 1200;

Opt.Method = 'matrix';
Opt.Symm = 'Ci';
Opt.nKnots = [30 0];

[x,y1] = salt(Sys1,Exp,Opt);
[x,y2] = salt(Sys2,Exp,Opt);

if opt.Display
  subplot(2,1,1)
  plot(x,y1,x,y2);
  legend('full Q','principal values');
  legend boxoff
  subplot(2,1,2);
  plot(x,y2-y1,'r');
  xlabel('frequency (MHz)');
  legend('difference');
  legend boxoff
end

ok = areequal(y1,y2,3e-3,'rel');
