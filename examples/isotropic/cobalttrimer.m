% Cobalt trimer
%==========================================================
% This demonstrates the difference between perturbation theory
% and the exact simulation of an isotropic solution spectrum of
% a Co trimer.
% See
% - Peake et al, Inorg. Chem. 18(4), 1000-1005 (1979)
% - P. Rieger, "Electron Spin Resonance", Royal Society of Chemistry,
%     2007, pp.48 (chapter 3.3 - Puzzling Line Shapes)

clear, clf

Sys.Nucs = 'Co';
Sys.A = 100;
Sys.n = 3;
Sys.lwpp = 0.8;

Exp.mwFreq = 9.5;

Opt.Method = 'perturb1';
[x,y1] = garlic(Sys,Exp,Opt);
Opt.Method = 'perturb2';
[x,y2] = garlic(Sys,Exp,Opt);
Opt.Method = 'exact';
[x,y0] = garlic(Sys,Exp,Opt);

plot(x,y1,x,y2,x,y0);
legend('1st order','2nd order','exact');
