% Naphthalene radical anion and cation
%==========================================================================
% Hyperfine values taken from Gerson book, Table 8.8, p. 234

clear, clf

Exp.mwFreq = 9.5;
Exp.Range = [337 342];
Exp.nPoints = 8192;

Sys.g = 2;
Sys.Nucs = '1H,1H';
Sys.n = [4 4];
Sys.lwpp = [0 ,0.01];  % Lorentzian lines

A_anion =  mt2mhz([-0.495 -0.183]);
A_cation = mt2mhz([-0.587 -0.167]);

Sys.A = A_anion;
[xa,ya] = garlic(Sys,Exp);
Sys.A = A_cation;
[xc,yc] = garlic(Sys,Exp);

plot(xa,ya,'b',xc,yc,'g');
legend('radical anion','radical cation');
xlabel('magnetic field [mT]');
set(gca,'YTick',[]);
