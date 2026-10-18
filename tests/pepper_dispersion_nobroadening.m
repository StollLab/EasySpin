function ok = test()

% Check that a dispersion signal without convolutional line width raises an
% error instead of returning the absorption spectrum, for powder (projection)
% and isotropic (template summation) spectra

Exp.mwFreq = 9.5;
Exp.Range = [315 345];
Exp.Harmonic = 0;
Exp.mwPhase = pi/2;

gs = {[2 2.05 2.1], 2.05};

for k = numel(gs):-1:1
  Sys.g = gs{k};
  try
    pepper(Sys,Exp);
    ok(k) = false;
  catch err
    ok(k) = contains(err.message,'Exp.mwPhase');
  end
end

end
