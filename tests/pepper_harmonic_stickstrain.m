function ok = test()

% Check the harmonic of a brute-force stick spectrum when only strains are
% given. The stick summation ignores the strains, so the default harmonic is
% 0, and an explicit first harmonic raises an error.

Sys.g = [2 2.05 2.1];
Sys.lwAnisoG = [50 100 150];
Exp.mwFreq = 9.5;
Exp.Range = [315 345];
Opt.BruteForce = true;

% Default harmonic
try
  pepper(Sys,Exp,Opt);
  ok(1) = true;
catch
  ok(1) = false;
end

% Explicit first harmonic
Exp.Harmonic = 1;
try
  pepper(Sys,Exp,Opt);
  ok(2) = false;
catch err
  ok(2) = contains(err.message,'No broadening given');
end

end
