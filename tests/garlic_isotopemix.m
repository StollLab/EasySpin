function ok = test()

% Isotope mixtures: stick spectrum integral should be independent of mixture

Exp.mwFreq = 9.668;
Exp.CenterSweep = [345 2];
Exp.Harmonic = 0;

Sys1.Nucs = 'B';
Sys1.A = 10;
[x,y1] = garlic(Sys1,Exp);

Sys1.Nucs = '11B';
Sys1.A = 10;
[x,y2] = garlic(Sys1,Exp);

ok = areequal(sum(y1),sum(y2),1e-6,'rel');
