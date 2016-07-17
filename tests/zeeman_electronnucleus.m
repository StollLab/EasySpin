function [err,data] = test(opt,olddata)

%======================================================
% Test 3: Electron and nucleus
%======================================================
Nuc = '14N';
Sys = struct('S',[1/2],'Nucs',Nuc,'g',[2 2.1 2.2],'A',[1 2 3]);
[Gx0,Gy0,Gz0] = zeeman(Sys,2);
f = -nucgval(Nuc)*nmagn/planck/1e9;

Gx1 = f*sop(Sys,'ex');
Gy1 = f*sop(Sys,'ey');
Gz1 = f*sop(Sys,'ez');

A =[Gx0 Gy0 Gz0];
B =[Gx1 Gy1 Gz1];
err = any(abs(A(:)-B(:))>1e-10);
data = [];
