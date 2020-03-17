function ok = test()

% Test zero-field energies of a and F terms against explicit equations
% (Bleaney/Trenam, Proc Roy Soc A, 223, 1-14 (1954))


Sys.S = 5/2;
D = 101;
a = 23;
F = 14.9;
Sys.D = D;
Sys.aF = [a F];
Sys.aFFrame = 3;

% Numerical energies
E1 = eig(zfield(Sys));
E1 = E1(1:2:6);

% Analytical energies
E2(1) = D/3 - (a-F)/2 - sqrt((18*D+a-F)^2+80*a^2)/6;
E2(2) = D/3 - (a-F)/2 + sqrt((18*D+a-F)^2+80*a^2)/6;
E2(3) = -2*D/3 + (a-F);
E2 = sort(E2.');

ok = areequal(E1,E2,1e-8,'abs');
