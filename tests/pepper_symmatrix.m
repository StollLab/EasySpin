function ok = test(opt)

% Symmetric matrices [xx yy zz xy xz yz] vs. full 3x3 matrices

mat2sym = @(M)[M(1,1) M(2,2) M(3,3) M(1,2) M(1,3) M(2,3)];

Rg = erot([0.2 0.5 -0.3]).';
RA = erot([-0.4 0.3 0.6]).';
gM = Rg*diag([2.0 2.1 2.2])*Rg.';
AM = RA*diag([20 40 100])*RA.';

Sys1.g = gM;
Sys1.Nucs = '1H';
Sys1.A = AM;
Sys1.lwpp = 1;

Sys2.g = mat2sym(gM);
Sys2.Nucs = '1H';
Sys2.A = mat2sym(AM);
Sys2.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.Range = [295 345];

[x,y1] = pepper(Sys1,Exp);
[x,y2] = pepper(Sys2,Exp);

if opt.Display
  plot(x,y1,x,y2,'r--');
  legend('full matrices','symmetric matrices');
  legend boxoff
end

ok = areequal(y1,y2,1e-10,'rel');
