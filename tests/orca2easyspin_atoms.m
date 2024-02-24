function ok = test(opt)

% Read binary ORCA property files

filename = 'orca/nx4me.prop';

Sys = orca2easyspin(filename);

NucsList = runprivate('nucstring2list',Sys.Nucs);

ok(1) = numel(Sys.data.NucId)==size(Sys.xyz,1);
ok(2) = numel(NucsList)==numel(Sys.NucsIdx);
