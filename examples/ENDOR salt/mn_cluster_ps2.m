% ENDOR of the tetranuclear Mn cluster in photosystem II (S2 state)
%----------------------------------------------------------------------------------------------------
clear, clf, clc

% g and A values from Peloquin et al, J.Am.Chem.Soc. 2000, 122, 10926
g = [1.97 1.99];
A = [311 270; 232 270; 200 250; 180 240];

% Experiment at W band - no need for mw frequency,
% since no orientation selection
Exp.Range = [0 250];
Exp.Field = 3380;
Exp.nPoints = 501;

PSII.g = g;
PSII.lwEndor = 3;

% (1) first-order perturbation treatment
PSII.Nucs = '55Mn,55Mn,55Mn,55Mn';
PSII.A = A;
Opt.Method = 'perturb1';
[freq,spec_p] = salt(PSII,Exp,Opt);
spec_p = spec_p/sum(spec_p);

% (2)matrix diagonalization, for each nucleus separately
PSII.Nucs = '55Mn';
Opt.Method = 'matrix';
for iMn=1:4
  PSII.A = A(iMn,:);
  [freq,spec_m(iMn,:)] = salt(PSII,Exp,Opt);
end
spec_m = sum(spec_m);
spec_m = spec_m/sum(spec_m);

plot(freq,spec_p,'r',freq,spec_m,'b');
xlabel('frequency (MHz)');
legend('Perturbation theory','matrix diagonalization');






