function ok = test()

% When Lanczos breaks down, chili falls back to Arnoldi
% (see https://github.com/StollLab/EasySpin/issues/178)

Sys.g = [2 2 2.1];
Sys.tcorr = 1e-9;

Exp.Field = 350;
Exp.mwRange = [9 11];
Exp.nPoints = 100;

Opt.LLMK = [8 0 2 0];

% Force Lanczos breakdown at the first step (|y.'*y| <= y'*y always holds)
Opt.Solver = 'L';
Opt.LanczosBreakdownTol = 2;
[~,spc] = chili(Sys,Exp,Opt);

Opt.Solver = '\';
[~,spcref] = chili(Sys,Exp,Opt);

ok = areequal(spc,spcref,1e-3,'rel');
