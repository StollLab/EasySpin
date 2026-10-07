function ok = test()

% Removal of spin-zero nuclei from per-nucleus and per-pair fields

% Nuclei 1H,12C,13C: 12C (I=0) is removed
% Pairs (1,2), (1,3), (2,3): only pair (1,3) remains
Sys.S = 1/2;
Sys.Nucs = '1H,12C,13C';
Sys.A = [1 1 1; 2 2 2; 3 3 3];
Sys.AFrame = [0.1 0 0; 0.2 0 0; 0.3 0 0];
Sys.sigma = [1 2 3; 4 5 6; 7 8 9]*1e-6 + 1;
Sys.sigmaFrame = [0 0.1 0; 0 0.2 0; 0 0.3 0];
Sys.nn = [10 11 12; 20 21 22; 30 31 32]*1e-3;
Sys.nnFrame = [0 0 0.1; 0 0 0.2; 0 0 0.3];

Q = runprivate('validatespinsys',Sys);
ok(1) = Q.nNuclei==2;
ok(2) = isequal(Q.A,Sys.A([1 3],:)) && isequal(Q.AFrame,Sys.AFrame([1 3],:));
ok(3) = isequal(Q.sigma,Sys.sigma([1 3],:)) && isequal(Q.sigmaFrame,Sys.sigmaFrame([1 3],:));
ok(4) = isequal(Q.nn,Sys.nn(2,:)) && isequal(Q.nnFrame,Sys.nnFrame(2,:));

% Full sigma and nn matrices
M = @(k) k*eye(3) + [0 1 2; 1 0 3; 2 3 0]*1e-3;
Sys = struct('S',1/2,'Nucs','1H,12C,13C','A',[1 1 1; 2 2 2; 3 3 3]);
Sys.sigma = [M(1); M(2); M(3)];
Sys.nn = [M(10); M(20); M(30)]*1e-3;
Q = runprivate('validatespinsys',Sys);
ok(5) = isequal(Q.sigma,[M(1); M(3)]) && isequal(size(Q.sigmaFrame),[2 3]);
ok(6) = isequal(Q.nn,M(20)*1e-3) && isequal(size(Q.nnFrame),[1 3]);

% Only one nucleus left: no nn couplings remain
Sys = struct('S',1/2,'Nucs','1H,12C','A',[1 1 1; 2 2 2],'nn',[1 2 3]*1e-3);
Q = runprivate('validatespinsys',Sys);
ok(7) = Q.nNuclei==1 && isempty(Q.nn) && isempty(Q.nnFrame);
