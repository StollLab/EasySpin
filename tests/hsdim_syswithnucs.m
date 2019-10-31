function err = test(opt,olddata)

% spin system with nuclei
%================================================================
Sys.S = 1/2;
Sys.g = [2 3 4];
Sys.Nucs = '14N';
Sys.A = [6 34 4];

N = hsdim(Sys);

Nref = 6;

err = ~(N==Nref);
