function ok = test()

% Legacy ORCA 3.0.3 files: check which fields are present in the spin
% system for calculations with various subsets of EPR properties, both for
% main output files (.oof) and binary property files (.prop).

folder = 'orca/v3.0.3/';

%   name              g      A      Q      D      Nucs  NucsIdx
calcs = {
  'dioxygen_g',       true,  false, false, false, '',   []
  'dioxygen_gD',      true,  false, false, true,  '',   []
  'hydroxyl_g',       true,  false, false, false, '',   []
  'hydroxyl_gA',      true,  true,  false, false, 'H',  2
  'hydroxyl_gAiso',   true,  true,  false, false, 'H',  2
  'hydroxyl_gQ',      true,  false, true,  false, 'H',  2
  'hydroxyl_Q',       false, false, true,  false, 'H',  2
  'hydroxyl_HO',      true,  true,  false, false, 'H',  1
  'hydroxyl_098_v303',true,  true,  true,  false, 'H',  2
  };

for k = size(calcs,1):-1:1
  [name,hasg,hasA,hasQ,hasD,Nucs,NucsIdx] = calcs{k,:};
  SysM = orca2easyspin([folder name '.oof']);
  SysP = orca2easyspin([folder name '.prop']);

  okk = true;
  for Sys = [{SysM},{SysP}]
    S = Sys{1};
    okk = okk && isfield(S,'g')==hasg && isfield(S,'A')==hasA;
    okk = okk && isfield(S,'Q')==hasQ && isfield(S,'D')==hasD;
    if isempty(Nucs)
      okk = okk && ~isfield(S,'Nucs');
    else
      okk = okk && ischar(S.Nucs) && strcmp(S.Nucs,Nucs) && isequal(S.NucsIdx,NucsIdx);
    end
  end

  % Spin is only stored in the main output file
  okk = okk && SysM.S==(1+startsWith(name,'dioxygen'))/2;

  ok(k) = okk;
end
