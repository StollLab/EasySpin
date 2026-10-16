function ok = test()

% Conversion of Sys.initState to a density matrix in the uncoupled basis

thresh = 1e-12;

% Electron spin 1/2 with one nucleus: dimensions and expansion
%-------------------------------------------------------------------------------
Sys = struct('S',1/2,'Nucs','1H','A',10);

rho0 = [0.4 0.1i 0 0.2; -0.1i 0.3 0 0; 0 0 0.2 0; 0.2 0 0 0.1];
Sys.initState = rho0;
rho = runprivate('initstate',runprivate('validatespinsys',Sys));
ok(1) = areequal(rho,rho0,thresh,'abs');

rhoe = [0.6 0.1; 0.1 0.4];
Sys.initState = {rhoe,'uncoupled'};
rho = runprivate('initstate',runprivate('validatespinsys',Sys));
ok(2) = areequal(rho,kron(rhoe,eye(2)/2),thresh,'abs');

Sys.initState = {[1 0],'uncoupled'};
rho = runprivate('initstate',runprivate('validatespinsys',Sys));
ok(3) = areequal(rho,diag([0.5 0.5 0 0]),thresh,'abs');

% Eigenbasis: lowest level of an S=1/2 at positive field is mS=-1/2
%-------------------------------------------------------------------------------
Sys = struct('S',1/2,'g',2);
Sys.initState = {[1 0],'eigen'};
Sys = runprivate('validatespinsys',Sys);
[rho,basis] = runprivate('initstate',Sys);
ok(4) = strcmp(basis,'eigen') && areequal(rho,diag([1 0]),thresh,'abs');
[rho,basis] = runprivate('initstate',Sys,ham(Sys,[0 0 350]));
ok(5) = strcmp(basis,'uncoupled') && areequal(rho,diag([0 1]),thresh,'abs');

% Triplet: zero-field basis vs xyz basis
% With D>0 and E>0, the zero-field levels are ordered Tz < Tx < Ty.
%-------------------------------------------------------------------------------
Sys = struct('S',1,'D',[1000 100]);
Sys.initState = {[0.2 0.5 0.3],'zerofield'};
rhoZF = runprivate('initstate',runprivate('validatespinsys',Sys));
Sys.initState = {[0.5 0.3 0.2],'xyz'};
rhoXYZ = runprivate('initstate',runprivate('validatespinsys',Sys));
ok(6) = areequal(rhoZF,rhoXYZ,thresh,'abs');

% Triplet with axial D: Tx and Ty are degenerate at zero field
%-------------------------------------------------------------------------------
Sys = struct('S',1,'D',1000);
Sys.initState = {[0.2 0.4 0.4],'zerofield'};
rhoZF = runprivate('initstate',runprivate('validatespinsys',Sys));
Sys.initState = {[0.4 0.4 0.2],'xyz'};
rhoXYZ = runprivate('initstate',runprivate('validatespinsys',Sys));
ok(7) = areequal(rhoZF,rhoXYZ,thresh,'abs');

Sys.initState = {[0.2 0.3 0.5],'zerofield'};
Sys = runprivate('validatespinsys',Sys);
try
  runprivate('initstate',Sys);
  ok(8) = false;
catch
  ok(8) = true;
end
