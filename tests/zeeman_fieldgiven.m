function [err,data] = test(opt,olddata)

%======================================================
% Test two forms of calls
%======================================================
Sys.S = 3/2;
Sys.g = diag([2 2 2])+rand/2;
B = ang2vec(rand*2*pi,rand*pi)*340;

[Gx0,Gy0,Gz0] = zeeman(Sys);
H0 = Gx0*B(1) + Gy0*B(2) + Gz0*B(3);
H1 = zeeman(Sys,B);

err = any(abs(H0(:)-H1(:))>1e-10);

data = [];
