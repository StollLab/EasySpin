function ok = test()

% Check that nonexistent and invalid files raise errors

cases = {
  'orca/v6.1.0/nonexistent.out',          'Cannot access'
  'orca/v6.1.0/nonexistent.property.txt', 'Cannot access'
  'orca/inputs/nitroxide.inp',            'not an ORCA output file'
  };

for k = size(cases,1):-1:1
  try
    orca2easyspin(cases{k,1});
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,cases{k,2});
  end
end
