function ok = test(opt)

% Read binary ORCA property files

filename = 'orca/nx4me';

Sys = orca2easyspin(filename);

Nucs = runprivate('nucstring2list',Sys.Nucs);

ok(1) = numel(Sys.Atoms)==size(Sys.xyz,1);
ok(2) = numel(Nucs)==numel(Sys.NucsIdx);
