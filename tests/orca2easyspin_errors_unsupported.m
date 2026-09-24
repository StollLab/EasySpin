function ok = test()

% Check that unsupported file formats raise errors: text-based property
% files from ORCA 4 and 5 (_property.txt)

versions = {'v4.2.1','v5.0.4'};
mols = {'aminoxyl','dioxygen','nitroxide','tripletformaldehyde'};

ok = true(numel(versions),numel(mols));
for v = 1:numel(versions)
  for m = 1:numel(mols)
    try
      orca2easyspin(['orca/' versions{v} '/' mols{m} '_property.txt']);
      ok(v,m) = false;
    catch err
      ok(v,m) = contains(err.message,'not supported');
    end
  end
end

ok = ok(:);
