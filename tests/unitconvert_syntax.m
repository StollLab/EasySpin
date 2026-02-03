function ok = test()

% Test whether unitconvert accepts all documented unit pairs

units = {'cm^-1', 'eV', 'J', 'K', 'MHz', 'GHz', 'THz', 'G', 'mT', 'T'};

nUnits = numel(units);

for u1 = 1:nUnits
  for u2 = 1:nUnits
    if u1==u2, continue; end
    unitpair = [units{u1} '->' units{u2}];
    [~] = unitconvert(rand,unitpair);
  end
end

ok = true;
