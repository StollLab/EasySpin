function ok = test()

% Non-equilibrium density matrix at negative fields. Reference: a rotation
% of the sample by pi around xL maps the field onto -zL, leaves the
% linearly polarized mw field along xL unchanged, and rotates the density
% matrix together with the molecule.

Sys.S = 1;
Sys.D = [300 50];
rho = [0.5 0.1+0.2i 0; 0.1-0.2i 0.3 0.05; 0 0.05 0.2];
Sys.initState = {rho,'uncoupled'};

Exp.mwFreq = 9.5;
Exp.SampleFrame = [0.2 0.7 0.1];

Exp.Range = [-500 -200];
[Pneg,Ineg] = resfields(Sys,Exp);

Exp.Range = [200 500];
Exp.SampleRotation = {[1;0;0],pi};
[Pref,Iref] = resfields(Sys,Exp);

ok = areequal(Pneg,-Pref,1e-6,'abs') && areequal(Ineg,Iref,1e-3,'rel');
