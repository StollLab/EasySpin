function [err,data] = test(opt,olddata)

% Test 2: sham vs zeeman/hfine
%======================================================
B = rand(1,3)*400;
Sys = struct('S',1/2,'g',[2 3 4],'Nucs','63Cu','A',[50 50 350]);
F = nquad(Sys) + hfine(Sys);
[Gx,Gy,Gz] = zeeman(Sys);
G = B(1)*Gx + B(2)*Gy + B(3)*Gz;
G = G/norm(B);
[F1,G1] = sham(Sys,B);

err = ~areequal(F1+norm(B)*G1,F+norm(B)*G);
data = [];
