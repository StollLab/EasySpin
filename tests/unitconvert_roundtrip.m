function ok = test()

% Test whether forward and backward conversions are inverse of each other.

units = {'cm^-1', 'eV', 'J', 'K', 'MHz', 'GHz', 'THz', 'G', 'mT', 'T'};

nUnits = numel(units);

ok = false(1,nUnits^2-nUnits);

p = 0;
for u1 = 1:nUnits
  for u2 = 1:nUnits
    if u1==u2, continue; end
    forward = [units{u1} '->' units{u2}];
    reverse = [units{u2} '->' units{u1}];
    v = unitconvert(unitconvert(1,forward),reverse);
    p = p+1;
    ok(p) = areequal(v,1,1e-10,'rel');
  end
end

end
