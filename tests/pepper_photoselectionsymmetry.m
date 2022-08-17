function ok = test()

% Test that triplet simulation with photoselection gives the same result
% for equivalent transition dipole moment orientations defined differently
% (correct result always obtained by setting GridSymmetry to C1)

% Define spin system
Triplet.S = 1;
Triplet.D = [-1000 150]; % MHz
Triplet.lwpp = 1; % mT
Triplet.initState = {[0.2 0.2 0.6],'zerofield'};  % non-equilibrium populations

Exp.mwFreq = 9.7;  % GHz
Exp.Range = [280 400];  % mT
Exp.Harmonic = 0;  % for transient EPR, turn off field modulation

Exp.lightScatter = 0;  % no isotropic contribution

Opt = struct();

% Transition dipole moment at 45 deg
Triplet.tdm = [0 45]*pi/180;  % orientation of tdm in molecular frame
    
Exp.lightBeam = 'parallel';
[~,spc_para1] = pepper(Triplet,Exp,Opt);

Exp.lightBeam = 'perpendicular';
[~,spc_perp1] = pepper(Triplet,Exp,Opt);
    
% Transition dipole moment at -45 deg
Triplet.tdm = -Triplet.tdm;
    
Exp.lightBeam = 'parallel';
[~,spc_para2] = pepper(Triplet,Exp,Opt);

Exp.lightBeam = 'perpendicular';
[~,spc_perp2] = pepper(Triplet,Exp,Opt);

ok(1) = areequal(spc_para1,spc_para2,1e-2,'rel');
ok(2) = areequal(spc_perp1,spc_perp2,1e-2,'rel');

end
