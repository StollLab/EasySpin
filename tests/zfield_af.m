function ok = test()

% Test a and F terms against explicit formulae and stev
% (Abragam/Bleaney p. 437)

a = 12;
F = 123;

S = 5/2;
Sys.S = S;

O40 = stev(S,[4 0]);
O44 = stev(S,[4 4]);

% (1) a term
%------------------------------------------------------------------
Sys.aF = [a 0];

H0 = a/120*(O40 + 5*O44);
H1 = zfield(Sys);

ok(1) = areequal(H1,H0,1e-6,'abs');

% (2) F term
%-----------------------------------------------------------------
Sys.aF = [0 F];

H0 = F/180*O40;
H1 = zfield(Sys);

ok(2) = areequal(H0,H1,1e-6,'abs');
