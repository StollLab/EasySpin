function ok = test()

% Test whether photoselect()'s integral over the third angle chi for unpolarized
% light is correct by comparing against explicit numerical integration.
% Rotations around chi do not affect resonance fields.

rng(11114); % seed rng for reproducibility

tdm = [rand*2*pi rand*pi]; % TDM orientation in mol frame
k = [rand*2*pi rand*pi]; % light propagation direction in lab frame

alpha = linspace(0,2*pi,3601); % all polarization angles

% Orientation of lab frame relative to mol frame; omitting chi to obtain
% integral over chi from photoselect()
ori = [rand*2*pi rand*pi];

% (1) Numerically integrate photoselection weights for all polarization angles
% alphe. Integral over chi is done by photoselect()
for ialpha = 1:numel(alpha)
  w(ialpha) = photoselect(tdm,ori,k,alpha(ialpha));
end
integratedWeight1 = mean(w);

% (2) Intregrals over alpha and chi are done by photoselect()
integratedWeight2 = photoselect(tdm,ori,k,NaN);

ok = areequal(integratedWeight1,integratedWeight2,1e-3,'rel');
