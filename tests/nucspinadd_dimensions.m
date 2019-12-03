function err = test(opt,olddata)

% Test 2: check matrix dimensions
%====================================================
A = [1 1 1];
AFrame = [2 3 4];
Q = [-1 0 1];
QFrame = [9 8 7];

Sys = struct('S',1/2,'g',[2 2 2.2]);

Sys = nucspinadd(Sys,'1H',A,AFrame,Q,QFrame);

if all(Sys.A(end,:)==A) && ...
   all(Sys.Q(end,:)==Q) && ...
   all(Sys.AFrame(end,:)==AFrame) && ...
   all(Sys.QFrame(end,:)==QFrame) && ...
   strcmp(Sys.Nucs(end-1:end),'1H')
  err = false;
else
  err = true;
end

