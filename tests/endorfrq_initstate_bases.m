function ok = test()

% Sys.initState: equivalent states given in different forms and bases give
% the same ENDOR frequencies and intensities. Populations enter only with
% orientation selection (Exp.mwFreq).

Exp.Field = 340;
Exp.Range = [0 60];
Exp.mwFreq = 9.9;
Exp.ExciteWidth = 500;
Opt.Threshold = 0;

% Population vector in the uncoupled basis vs diagonal density matrix,
% electron-only state vs full state
Sys = struct('S',1/2,'g',[2 2.1 2.2],'Nucs','1H','A',[10 20 30]);
pop = [0.7 0.3];
Sys.initState = {pop,'uncoupled'};
[nu1,I1] = endorfrq(Sys,Exp,Opt);
Sys.initState = diag(kron(pop,[1 1]/2));
[nu2,I2] = endorfrq(Sys,Exp,Opt);
ok(1) = areequal(nu1,nu2,1e-10,'rel') && areequal(I1,I2,1e-10,'rel');

% Eigenbasis vs uncoupled basis, with the field along zMol: the two lowest
% levels belong to mS = -1/2. The small hyperfine mixing of the mS = +-1/2
% states changes the populations only in second order.
Sys = struct('S',1/2,'g',2,'Nucs','1H','A',10);
Exp.mwFreq = 9.52;
Sys.initState = {[0.4 0.4 0.1 0.1],'eigen'};
[nu1,I1] = endorfrq(Sys,Exp,Opt);
Sys.initState = {[0.2 0.8],'uncoupled'};
[nu2,I2] = endorfrq(Sys,Exp,Opt);
ok(2) = areequal(nu1,nu2,1e-10,'rel') && areequal(I1,I2,1e-5,'rel');

% Triplet: zero-field basis vs xyz basis. With D>0 and E>0, the zero-field
% levels are ordered Tz < Tx < Ty. The weak hyperfine coupling mixes the
% zero-field states only slightly.
Sys = struct('S',1,'D',[900 160],'Nucs','1H','A',[2 3 4]);
Exp.mwFreq = 9.5;
Sys.initState = {[0.1 0.5 0.4],'zerofield'};
[nu1,I1] = endorfrq(Sys,Exp,Opt);
Sys.initState = {[0.5 0.4 0.1],'xyz'};
[nu2,I2] = endorfrq(Sys,Exp,Opt);
ok(3) = areequal(nu1,nu2,1e-10,'rel') && areequal(I1,I2,1e-3,'rel');

% Scaling: intensities are proportional to the state, and only population
% differences matter, so a traceless state works too
Sys = struct('S',1/2,'g',[2 2.1 2.2],'Nucs','1H','A',[10 20 30]);
Exp.mwFreq = 9.9;
Sys.initState = {[0.7 0.3],'uncoupled'};
[~,I1] = endorfrq(Sys,Exp,Opt);
Sys.initState = {[1.4 0.6],'uncoupled'};
[~,I2] = endorfrq(Sys,Exp,Opt);
Sys.initState = {[0.2 -0.2],'uncoupled'};
[~,I3] = endorfrq(Sys,Exp,Opt);
ok(4) = areequal(I2,2*I1,1e-10,'rel') && areequal(I3,I1,1e-10,'rel');
