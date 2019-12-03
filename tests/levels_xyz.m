function err = test(opt,olddata)

% Shorthands 'x', 'y', 'z' for orientations

Sys = struct('S',1,'D',10e3);
B = linspace(0,100,201);
Ex = levels(Sys,'x',B);
Ey = levels(Sys,'y',B);
Ez = levels(Sys,'z',B);
Exy = levels(Sys,'xy',B);
Eyz = levels(Sys,'yz',B);
Exz = levels(Sys,'xz',B);
Exyz = levels(Sys,'xyz',B);

err = false;
