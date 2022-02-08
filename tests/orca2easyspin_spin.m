function ok = test()

% Check for correct spin multiplicity
%-------------------------------------------------
% Orca <=3.0.3 does not save spin multiplicity in the .prop file.
% Orca 4.0.0 does.

ok = true;

% test doublet and triplet
FileName{1} = 'orca/hydroxyl_098_orca400.oof';
FileName{2} = 'orca/dioxygen_gD_orca400.oof';
spin = [1/2, 1];

for f = 1:numel(FileName)
  Sys = orca2easyspin(FileName{f});
  ok = ok && Sys.S==spin(f);
end
