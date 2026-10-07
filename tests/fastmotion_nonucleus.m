function ok = test()

% Linewdiths of hypothetical system without nuclear spins
%======================================================

Sys.g = [2.0 2.1 2.2];
Field = 350;
tcorr = 1e-11;

lw = fastmotion(Sys,Field,tcorr);

lw_correct = 0.41832353;

ok(1) = areequal(lw,lw_correct,1e-5,'rel');

% isotropic g without nuclei gives an informative error
Sys.g = 2;
try
  fastmotion(Sys,Field,tcorr);
  ok(2) = false;
catch e
  ok(2) = contains(e.message,'anisotropic');
end
