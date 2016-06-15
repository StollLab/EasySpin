function [err,data] = test(opt,olddata)

%======================================================
% Test 2: Two electrons
%======================================================
Sys = struct('S',[1/2 1],'g',[2 2.1 2.2; 3 4 5],'ee',[1 2 3]);
[Gx0,Gy0,Gz0] = zeeman(Sys,2);
f = Sys.g(2,:)*bmagn/planck/1e9;

Gx1 = f(1)*sop(Sys,'ex');
Gy1 = f(2)*sop(Sys,'ey');
Gz1 = f(3)*sop(Sys,'ez');

A =[Gx0 Gy0 Gz0];
B =[Gx1 Gy1 Gz1];
err = any(abs(A(:)-B(:))>1e-10);
data = [];
