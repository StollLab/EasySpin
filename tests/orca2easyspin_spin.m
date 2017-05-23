function [err,data] = test(opt,olddata)

% Check for correct spin multiplicity
%-------------------------------------------------
% Orca <=3.0.3 does not save spin multiplicity in the .prop file.
% Orca 4.0.0 does.

err = false;

% test doublet and triplet
FileName{1} = 'orca/hydroxyl_098_orca400';
FileName{2} = 'orca/dioxygen_gD_orca400';
spin = [1/2, 1];

for f = 1:numel(FileName)
  Sys = orca2easyspin(FileName{f});
  err = err || Sys.S~=spin(f);
end

data = [];
