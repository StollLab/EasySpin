function ok = test()

% Check against explicitly expressions

units = {'cm^-1', 'eV', 'K', 'MHz', 'GHz', 'THz', 'G', 'mT', 'T'};

value_in = 1;

n = numel(units);
ok = false(1,n^2-n);
p = 0;
for u1 = 1:n
  for u2 = 1:n
    if u1==u2, continue; end

    % Convert value from source unit to J
    switch units{u1}
      case 'cm^-1'
        v_J = value_in*100*planck*clight;
      case 'eV'
        v_J = value_in*evolt;
      case 'J'
        v_J = value_in;
      case 'K'
        v_J = value_in*boltzm;
      case 'MHz'
        v_J = value_in*1e6*planck;
      case 'GHz'
        v_J = value_in*1e9*planck;
      case 'THz'
        v_J = value_in*1e12*planck;
      case 'G'
        v_J = value_in*1e-4*bmagn*gfree;
      case 'mT'
        v_J = value_in*1e-3*bmagn*gfree;
      case 'T'
        v_J = value_in*bmagn*gfree;
    end

    % Convert value from J to target unit
    switch units{u2}
      case 'cm^-1'
        value_ref = v_J/planck/clight/100;
      case 'eV'
        value_ref = v_J/evolt;
      case 'J'
        value_ref = v_J;
      case 'K'
        value_ref = v_J/boltzm;
      case 'MHz'
        value_ref = v_J/planck/1e6;
      case 'GHz'
        value_ref = v_J/planck/1e9;
      case 'THz'
        value_ref = v_J/planck/1e12;
      case 'G'
        value_ref = v_J/bmagn/gfree/1e-4;
      case 'mT'
        value_ref = v_J/bmagn/gfree/1e-3;
      case 'T'
        value_ref = v_J/bmagn/gfree;
    end

    % Calculate value using unitconvert
    fromto = [units{u1} '->' units{u2}];
    value_out = unitconvert(value_in,fromto);

    % Compare
    p = p+1;
    ok(p) = areequal(value_out,value_ref,1e-8,'rel');
  end

end
