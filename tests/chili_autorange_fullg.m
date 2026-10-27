function ok = test()

% Automatic field range with a full g matrix is the same as with its
% principal values

Sys.g = [2.002 2.006 2.009];
Sys.tcorr = 1e-9;
Exp.mwFreq = 9.5;
Exp.nPoints = 200;

[B0,~] = chili(Sys,Exp);

R = erot(0.3,0.7,-0.2);
Sys.g = R*diag(Sys.g)*R.';
[B,~] = chili(Sys,Exp);

ok = areequal(B,B0,1e-10,'rel');

end
