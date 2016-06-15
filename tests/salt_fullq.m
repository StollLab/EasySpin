function [err,data] = test(opt,olddata)

% full Q tensor

Qpv = [-1 -1 2]*0.1;
QFrame = [10 30 10]*pi/180;
RQ = erot(QFrame);

Sys1.Nucs = '14N';
Sys1.A = 1;
Sys1.Q = RQ.'*diag(Qpv)*RQ;
Sys1.QFrame = [0 0 0];

Sys2.Nucs = '14N';
Sys2.A = 1;
Sys2.Q = Qpv;
Sys2.QFrame = QFrame;

Sys1.lwEndor = 0.1;
Sys2.lwEndor = 0.1;

Exp.Range = [0 7];
Exp.Field = 1200;

Opt.Method = 'matrix';
Opt.Symm = 'Ci';

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

err = ~areequal(y1/max(y1),y2/max(y2),1e-3);

data = [];
