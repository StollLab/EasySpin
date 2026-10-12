function ok = test()

% Left-to-right (Lentz) and right-to-left continued-fraction evaluation
% must give the same spectrum on the full abscissa

Sys.g = [2.008 2.006 2.003];
Sys.tcorr = 1e-9;

Exp.Range = [336 340];
Exp.mwFreq = 9.4788;
Exp.nPoints = 1000;

Opt.LLMK = [8 0 2 0];

Opt.Lentz = true;
[~,spc1] = chili(Sys,Exp,Opt);
Opt.Lentz = false;
[~,spc2] = chili(Sys,Exp,Opt);

ok = numel(spc2)==Exp.nPoints && areequal(spc1,spc2,1e-4,'rel');
