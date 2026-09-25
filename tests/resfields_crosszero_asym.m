function ok = test()

% Asymmetric field range across zero: mirrored resonances partly in range
% are returned with all their positions (no NaN)

Sys.g = [2 2.1 2.2];
Exp.mwFreq = 9.5;
Exp.SampleFrame = [0 0 0; 0 pi/2 0; pi/2 pi/2 0];  % 339, 323, 309 mT

Exp.Range = [0 400];
P0 = resfields(Sys,Exp);
Exp.Range = [-320 400];
P = resfields(Sys,Exp);

ok = areequal(P,[P0;-P0],1e-10,'abs') && ~any(isnan(P(:)));
