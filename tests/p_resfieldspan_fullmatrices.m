function ok = test()

% Rotated full g and A matrices give the same field span as their principal
% values, also with several nuclei

mwFreq = 9.5;  % GHz

Sys.g = [2.002 2.006 2.009];
Sys.Nucs = '14N,1H';
Sys.A = [15 15 95; -3 5 12];  % MHz
[minB0,maxB0] = runprivate('p_resfieldspan',runprivate('validatespinsys',Sys),mwFreq);

R1 = erot(0.3,0.7,-0.2);
R2 = erot(-1.1,0.4,0.9);
Sys.g = R1*diag(Sys.g)*R1.';
Sys.A = [R1*diag(Sys.A(1,:))*R1.'; R2*diag(Sys.A(2,:))*R2.'];
[minB,maxB] = runprivate('p_resfieldspan',runprivate('validatespinsys',Sys),mwFreq);

% Reference: first-order estimate from principal values
hf = (1*95 + 1/2*12)*1e6;  % Hz
minBref = planck*(mwFreq*1e9 - hf)/bmagn/2.009/1e-3;  % mT
maxBref = planck*(mwFreq*1e9 + hf)/bmagn/2.002/1e-3;  % mT

tol = 1e-10;
ok(1) = areequal(minB,minB0,tol,'rel');
ok(2) = areequal(maxB,maxB0,tol,'rel');
ok(3) = areequal([minB0 maxB0],[minBref maxBref],tol,'rel');

end
