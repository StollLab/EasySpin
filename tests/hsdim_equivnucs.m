function err = test(opt,olddata)

% spin system with equivalent nuclei
%================================================================
Sys1.S = 1/2;
Sys1.Nucs = '55Mn';
Sys1.n = 2;
Sys1.A = [5 6 7];

Sys2.S = 1/2;
Sys2.Nucs = '55Mn,55Mn';
Sys2.A = [5 6 7; 5 6 7];

N1 = hsdim(Sys1);
N2 = hsdim(Sys2);

err = ~areequal(N1,N2);
