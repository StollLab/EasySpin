function ok = test()

% Test whether photoselect()'s integral over the third angle chi for polarized
% light is correct by comparing against explicit numerical integration.
% Rotations around chi do not affect resonance fields.

rng(11114); % seed rng for reproducibility

tdm = [rand*2*pi rand*pi]; % TDM orientation in mol frame
k = [rand*2*pi rand*pi]; % light propagation direction in lab frame
alpha = rand*2*pi; % polarization angle

% Orientation of lab frame relative to mol frame
phi = rand*2*pi;
theta = rand*pi;
chi = linspace(0,2*pi,3601); % third angle, does not affect resonance field

% (1) Calculate integral numerically
for ichi = 1:numel(chi)
  ori(ichi,:) = [phi theta chi(ichi)];
end
w = photoselect(tdm,ori,k,alpha);
integratedWeight1 = sum(w)/numel(chi);

% (2) Calculate integral using photoselect()
integratedWeight2 = photoselect(tdm,[phi theta],k,alpha);

ok = areequal(integratedWeight1,integratedWeight2,1e-3,'rel');
