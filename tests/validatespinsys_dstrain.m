function ok = test()

% D strain input: 1 or 2 columns accepted, more rejected

Sys = struct('S',1,'D',[300 50]);

Sys.DStrain = 30;
[Sys_,err] = runprivate('validatespinsys',Sys);
ok(1) = isempty(err) && isequal(Sys_.DStrain,[30 0]);

Sys.DStrain = [30 10];
[Sys_,err] = runprivate('validatespinsys',Sys);
ok(2) = isempty(err) && isequal(Sys_.DStrain,[30 10]);

Sys.DStrain = [30 10 0.5];
[~,err] = runprivate('validatespinsys',Sys);
ok(3) = contains(err,'DStrainCorr');
