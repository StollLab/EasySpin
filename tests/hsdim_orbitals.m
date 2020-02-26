function ok = test()

% system with orbital angular momenta

Sys.S = 3/2;
Sys.Nucs = '1H';
Sys.L = 2;
Sys.A = [5 12 45];
Sys.soc = 1000;

N = hsdim(Sys);

S = Sys.S;
I = nucspin(Sys.Nucs);
L = Sys.L;

Nref = prod(2*[S I L]+1);

ok = (N==Nref);
