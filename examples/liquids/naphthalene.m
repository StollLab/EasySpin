% Naphthalene radical anion and cation
%==========================================================================
% Hyperfine values taken from Gerson book, Table 8.8, p. 234

clear, clf

Exp.mwFreq = 9.5;       % GHz
Exp.Range = [337 342];  % mT
Exp.nPoints = 8192;

Naph.g = 2;
Naph.Nucs = '1H,1H';
Naph.n = [4 4];
Naph.lwpp = [0, 0.01];  % Lorentzian linshape, mT

A_anion =  unitconvert([-0.495 -0.183],'mT->MHz');
A_cation = unitconvert([-0.587 -0.167],'mT->MHz');

Naph.A = A_anion;
[B,spc_anion] = garlic(Naph,Exp);
Naph.A = A_cation;
[B,spc_cation] = garlic(Naph,Exp);

plot(B,spc_anion,B,spc_cation);
legend('radical anion','radical cation');
xlabel('magnetic field [mT]');
set(gca,'YTick',[]);
