function [angles_M2L,nOrientations,nSites,averageOverChi] = p_crystalorientations(Exp,Opt)

% Function that takes care of crystal frame setups for single-crystal simulation
%    and a few related things.
% Used in resfields*, resfreqs*, endorfrq*, curry, saffron, resfields_eig

%{
Inputs:
  Opt.crystalSample       true/false (defined if coming from pepper/salt/saffron/curry)
  Exp.CrystalSymmetry     space group, only used for crystal simulations
  Exp.MolFrame            set by pepper/salt/saffron for powders, by user for cystals
  Opt.Sites               list of sites to include, e.g. [1 3] in a 4-site group
  Opt.R_L2S               transformation from lab to sample frame

Outputs:
  angles_M2L       list of orientations Euler angles, one per row
                     [phi,theta,chi] for mol-to-lab frame transformation
  nOrientations    total number of orientations
  nSites           number of symmetry-related sites
  averageOverChi   whether to compute average over chi angle
%}

% MolFrame contains an Nx3 array of Euler angles (in radians)
% specifying the relative orientation between sample and molecular frame.
if size(Exp.MolFrame,2)~=3
  error('Exp.MolFrame requires three numbers per row (three Euler angles).')
end

pepperCall = isfield(Opt,'rotatedSample');

% Process sample orientation if not already done by
% pepper/salt/saffron/curry
if ~pepperCall
  [Opt.R_L2S,Opt.rotatedSample] = p_samplerotation(Exp);
  Opt.crystalSample = true;
end

if Opt.crystalSample
  
  % Crystals
  %-----------------------------------------------------------------------------
  % Get site-to-site transformations for given crystal space group
  % (matrices for active rotations; defined in sample/crystal frame)
  if isempty(Exp.CrystalSymmetry)
    Exp.CrystalSymmetry = 'P1';
    logmsg(1,'  no crystal symmetry given: assuming space group P1 (#1)');
  end
  Rsite_C = sitetransforms(Exp.CrystalSymmetry);
  nSites  = numel(Rsite_C);
  logmsg(1,'  crystal symmetry: %d magnetically distinct sites in unit cell',nSites);
  
  % Select sites if user wants it
  if isfield(Opt,'Sites') && ~isempty(Opt.Sites)
    if any(Opt.Sites<1) || any(Opt.Sites>nSites)
      error('For the given crystal symmetry, Opt.Sites must contain numbers between 1 and %d.',nSites);
    end
    Rsite_C = Rsite_C(Opt.Sites);
    nSites = numel(Rsite_C);
  end
  
  % Sample/crystal-to-molecular frame transformation, R_S2M
  if ~isempty(Exp.MolFrame)
    logmsg(1,'  molecular frame given in Exp.MolFrame');
    if ~numel(Exp.MolFrame)==3
      error('Exp.MolFrame has wrong size.');
    end
    R_S2M = erot(Exp.MolFrame);
  else
    R_S2M = eye(3);
  end
  R_M2S = R_S2M.';

  % Generate list of mol frame orientations, represented in sample/crystal frame
  xyzM_M = eye(3);           % xM, yM, zM in mol frame representation
  xyzM0_S = R_M2S*xyzM_M;   % transformation to crystal/sample frame
  for iSite = nSites:-1:1
    xyzM_S{iSite} = Rsite_C{iSite}*xyzM0_S;   % active rotation in sample/crystal frame
  end

  % Generate list of lab frame orientations, represented in the
  % various site molecular frames
  % (M = molecular frame of reference site)
  % (Mi = molecular frame of site # i)
  nSampleOrientations = numel(Opt.R_L2S);
  angles_M2L = zeros(nSampleOrientations*nSites,3);
  idx = 1;
  for iSampleOri = 1:nSampleOrientations
    R_S2L_ = Opt.R_L2S{iSampleOri}.';
    for iSite = 1:nSites
      xyzMi_L = R_S2L_*xyzM_S{iSite};  % transformation to lab frame
      angles_M2L(idx,:) = eulang(xyzMi_L,true);  % Euler angles for Mi -> L frame
      idx = idx + 1;
    end
  end
  nOrientations = nSampleOrientations*nSites;
  
  averageOverChi = false;
  
else
  
  % Disordered and partially ordered systems (require grid integration)
  %-----------------------------------------------------------------------------
  % Grid orientations are supplied by pepper etc in Exp.MolFrame.
  
  % Exp.MolFrame transforms from sample(=lab) to molecular frame, but
  % angles_M2L = [phi, theta, chi] from mol to lab frame.
  % (Exp.SampleFrame is ignored for a powder)
  angles_M2L = -fliplr(Exp.MolFrame);
  nOrientations = size(angles_M2L,1);
  
  % For powder simulations, always average over the chi angle.
  averageOverChi = true;
  
  % Apply sample rotations(s)
  % Exp.SampleOrientation is disregarded for the purpose of grid construction.
  % (but needs to be taken into account if Exp.Ordering is anisotropic.)
  
  % Exp.CrystalSymmetry and Exp.MolFrame are disregarded, since they have
  % no effect on the powder spectrum.
  nSites = 1;
  if ~isempty(Exp.CrystalSymmetry)
    logmsg(1,'  powder simulation; disregarding Exp.CrystalSymmetry');
  end
  if ~isempty(Exp.MolFrame)
    logmsg(1,'  powder simulation; disregarding Exp.MolFrame');
  end

end
