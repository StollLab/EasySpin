function [err,data] = test(opt,olddata)

% Test 2: check matrix dimensions
%====================================================
A = [1 1 1];
AFrame = [0 0 0];
Q = [-1 0 1];
QFrame = [0 0 0];
Sys = struct('S',1/2,'g',[2 2 2.2]);
Sys = nucspinadd(Sys,'1H',A,AFrame,Q,QFrame);

if (Sys.A(end,:)==A) & ...
   (Sys.Q(end,:)==Q) & ...
   (Sys.AFrame(end,:)==AFrame) & ...
   (Sys.QFrame(end,:)==QFrame) & ...
   strcmp(Sys.Nucs(end-1:end),'1H')
  err = 0;
else
  err = 1;
end
data = [];
