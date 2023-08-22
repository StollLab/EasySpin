function ok = test()

% Test whether the space group information provided in Exp.CrystalSymmetry
% yield the expected resonance fields for a crystal by comparing to a
% explicit calculation that manually transforms all the frames involved.

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [300 400];  % mT

Sys.g = [2 2.1 2.2];

Sys.gFrame = [10 20 30]*pi/180;
Exp.MolFrame = [15 25 35]*pi/180;
Exp.SampleFrame = [38 48 67]*pi/180;

% Transformation matrices
R_L2S = erot(Exp.SampleFrame);
R_S2M = erot(Exp.MolFrame);
R_M2g = erot(Sys.gFrame);

nB0_L = [0; 0; 1];    % B0 is along zLab
nB0_S = R_L2S*nB0_L;
g_g = diag(Sys.g);

% Test one space group from each of the 11 Laue classes
spaceGroups = [2 15 74 88 142 148 167 176 194 206 230];
for sg = 1:numel(spaceGroups)
  Exp.CrystalSymmetry = spaceGroups(sg);
  Bres = resfields(Sys,Exp);

  % Get list of site rotation matrices for the space groups
  Rsite = runprivate('sitetransforms',spaceGroups(sg));
  nSites = numel(Rsite);

  % Determine nB0 vector in g frame for each site
  nB0_g = zeros(3,nSites);
  for s = 1:nSites
    R_S2Mi = R_S2M*Rsite{s}.';
    nB0_Mi = R_S2Mi*nB0_S;
    nB0_g(:,s) = R_M2g*nB0_Mi;
  end

  % Calculate effective g values and resonance fields
  geff = vecnorm(nB0_g.'*g_g,2,2);
  Bres_ref = planck*Exp.mwFreq*1e9/bmagn./geff/1e-3;  % mT

  ok(sg) = all(abs(Bres-Bres_ref)<1e-4);
end
