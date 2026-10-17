function ok = test()

% Sys.initState: equivalent states given in different forms and bases give
% the same resonance frequencies and intensities

Exp.Field = 340;
Opt.Threshold = 0;

% Population vector in the uncoupled basis vs diagonal density matrix,
% electron-only state vs full state
Sys = struct('S',1/2,'g',[2 2.1 2.2],'Nucs','1H','A',[10 20 30]);
pop = [0.7 0.3];
Sys.initState = {pop,'uncoupled'};
[nu1,I1] = resfreqs_matrix(Sys,Exp,Opt);
Sys.initState = diag(kron(pop,[1 1]/2));
[nu2,I2] = resfreqs_matrix(Sys,Exp,Opt);
ok(1) = areequal(nu1,nu2,1e-10,'rel') && areequal(I1,I2,1e-10,'rel');

% Eigenbasis vs uncoupled basis, with the field along zMol: the lowest
% level is mS = -1/2
Sys = struct('S',1/2,'g',2);
Sys.initState = {[0.8 0.2],'eigen'};
[nu1,I1] = resfreqs_matrix(Sys,Exp,Opt);
Sys.initState = {[0.2 0.8],'uncoupled'};
[nu2,I2] = resfreqs_matrix(Sys,Exp,Opt);
ok(2) = areequal(nu1,nu2,1e-10,'rel') && areequal(I1,I2,1e-10,'rel');

% Higher-order spin Hamiltonian: zero-field basis vs density matrix
% constructed from the zero-field eigenstates
Sys = struct('S',1,'Ham022',[100 0 300 0 50]);
Sys.Ham110 = -sqrt(3)*2*bmagn/(planck*1e9);  % g = 2
pop = [0.1 0.5 0.4];
Sys.initState = {pop,'zerofield'};
[nu1,I1] = resfreqs_matrix(Sys,Exp,Opt);
[V,E] = eig(ham(Sys,zeros(1,3)));
[~,idx] = sort(real(diag(E)));
V = V(:,idx);
Sys.initState = V*diag(pop)*V';
[nu2,I2] = resfreqs_matrix(Sys,Exp,Opt);
ok(3) = areequal(nu1,nu2,1e-10,'rel') && areequal(I1,I2,1e-8,'rel');
