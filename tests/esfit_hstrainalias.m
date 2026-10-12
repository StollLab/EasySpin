function ok = test()

% esfit renames Sys.HStrain (older name) to Sys.lwAnisoG once before fitting,
% so the renaming warning is shown only once

Sys.g = [2 2.1];
Sys.lwAnisoG = 15;

Exp.Field = 350;
Exp.mwRange = [9.5 10.5];

spc = pepper(Sys,Exp);

FitOpt.Verbosity = 0;
FitOpt.GridSize = 3;
FitOpt.RandomizeGrid = false;
FitOpt.Method = 'grid fcn';

Sys0 = rmfield(Sys,'lwAnisoG');
Sys0.HStrain = 10;
Vary.HStrain = 20;  % lower bound is clipped to zero

% with Vary
out = evalc('result = esfit(spc,@pepper,{Sys0,Exp},{Vary},FitOpt);');
nWarnings = numel(strfind(out,'Sys.HStrain has been renamed'));
ok(1) = nWarnings==1;
ok(2) = isfield(result.argsfit{1},'lwAnisoG') && ~isfield(result.argsfit{1},'HStrain');
ok(3) = result.rmsd/max(result.fit) < 1e-10;

% with lower and upper bounds
lb = Sys0; lb.HStrain = 0;
ub = Sys0; ub.HStrain = 30;
out = evalc('result = esfit(spc,@pepper,{Sys0,Exp},{lb},{ub},FitOpt);');
nWarnings = numel(strfind(out,'Sys.HStrain has been renamed'));
ok(4) = nWarnings==1;
ok(5) = result.rmsd/max(result.fit) < 1e-10;
