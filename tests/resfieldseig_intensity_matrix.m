function ok = test()

%=======================================================
% resfields_eig should give the same integral intensity as matrix diagonalization
% (tested for isotropic powder and single crystal)
%=======================================================
clear Sys Exp

Sys.g = 2.1;
Sys.Nucs = '1H';
Sys.A = 100;
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.Range = [300 360];
Exp.Harmonic = 0;

% isotropic powder
Opt.Method = 'matrix';
[x,y1] = pepper(Sys,Exp,Opt);
Opt.Method = 'eig';
[x,y2] = pepper(Sys,Exp,Opt);
dx = x(2)-x(1);

integral_matrix = sum(y1)*dx;
integral_resfields_eig = sum(y2)*dx;

ok(1) = areequal(integral_matrix,integral_resfields_eig,0.001,'abs');

% single crystal
Exp.CrystalOrientation = rand(1,3)*pi;
Opt.Method = 'matrix';
[x,y1] = pepper(Sys,Exp,Opt);
Opt.Method = 'eig';
[x,y2] = pepper(Sys,Exp,Opt);
dx = x(2)-x(1);

integral_matrix = sum(y1)*dx;
integral_resfields_eig = sum(y2)*dx;

ok(2) = areequal(integral_matrix,integral_resfields_eig,0.001,'abs');
