function ok = test(opt)

%=======================================================
% resfields_eig should give the same integral intensity as resfields
% (tested using pepper simulation of isotropic powder and single crystal)
%=======================================================

Sys.g = 2.1;
Sys.Nucs = '1H';
Sys.A = 100;
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.Range = [300 360];
Exp.Harmonic = 0;

% isotropic powder
Opt.Method = 'matrix';
[B,spc1] = pepper(Sys,Exp,Opt);
Opt.Method = 'eig';
[B,spc2] = pepper(Sys,Exp,Opt);
dx = B(2)-B(1);

integral_matrix = sum(spc1)*dx;
integral_resfields_eig = sum(spc2)*dx;

ok(1) = areequal(integral_matrix,integral_resfields_eig,0.001,'abs');

if opt.Display
  subplot(2,1,1)
  plot(B,spc1,B,spc2)
  title('isotropic powder')
end

% single crystal
Exp.MolFrame = [30 10 -20]*pi/180;
Opt.Method = 'matrix';
[B,spc1] = pepper(Sys,Exp,Opt);
Opt.Method = 'eig';
[B,spc2] = pepper(Sys,Exp,Opt);
dx = B(2)-B(1);

integral_matrix = sum(spc1)*dx;
integral_resfields_eig = sum(spc2)*dx;

ok(2) = areequal(integral_matrix,integral_resfields_eig,0.001,'abs');

if opt.Display
  subplot(2,1,2)
  plot(B,spc1,B,spc2)
  title('crystal')
end
