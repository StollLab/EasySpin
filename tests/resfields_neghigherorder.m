function ok = test()

% Negative fields with higher-order Zeeman terms: resonance fields agree
% with direct diagonalization at the negative field

Sys.S = 3/2;
Sys.g = 2;
Sys.D = 500;
Sys.Ham314 = 1e-6*[1 0 2 0 3 0 1 0 2];
Exp.mwFreq = 9.5;
Exp.Range = [-500 -100];

P = resfields(Sys,Exp);

ok = ~isempty(P);
for k = 1:numel(P)
  E = eig(ham(Sys,[0 0 P(k)]));
  dE = abs(E-E.');
  ok = ok && min(abs(dE(:)-Exp.mwFreq*1e3))<1e-2;  % MHz
end
