function ok = test()

% Check whether transition intensity follows crystal rotation

Sys.S = 1;
Sys.D = 30e3*2*[1 0.1];
Opt.Transitions = [1 2; 1 3; 2 3];

% Rotation of crystal around z, compensated by B1 rotation around z
Exp.CrystalOrientation = [0 0 0];
Exp.mwMode = {0 0};
[~,i0] = resfreqs_matrix(Sys,Exp,Opt);

a = 0.234543*pi;
Exp.CrystalOrientation = [a 0 0];
Exp.mwMode = {0 -a};
[~,i1] = resfreqs_matrix(Sys,Exp,Opt);

ok = areequal(i0,i1,1e-8,'abs');
