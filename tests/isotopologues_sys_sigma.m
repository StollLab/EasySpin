function ok = test()

% Isotopologues with sigma and sigmaFrame: rows of spin-0 isotopes are removed

Sys.Nucs = 'C,1H';
Sys.A = [5 10];
Sys.sigma = [1.001 1.002];
Sys.sigmaFrame = [0.1 0.2 0.3; 0.4 0.5 0.6];

Iso = isotopologues(Sys);
i12 = strcmp({Iso.Nucs},'1H');
i13 = strcmp({Iso.Nucs},'13C,1H');

ok(1) = isequal(Iso(i12).sigma,1.002) && isequal(Iso(i12).sigmaFrame,[0.4 0.5 0.6]);
ok(2) = isequal(Iso(i13).sigma,[1.001; 1.002]) && isequal(Iso(i13).sigmaFrame,Sys.sigmaFrame);

% Full 3x3 sigma matrices
Sys = struct('Nucs','1H,C,1H','A',[1 2 3],'sigma',kron([1;2;3],eye(3)));
Iso = isotopologues(Sys);
i12 = strcmp({Iso.Nucs},'1H,1H');
ok(3) = isequal(Iso(i12).sigma,kron([1;3],eye(3)));
