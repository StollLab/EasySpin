function [angles_M2L,nOrientations,nSites,averageOverChi] = p_crystalorientations(Exp,Opt)

% Function that takes care of crystal frame setups for single-crystal simulation
%    and a few related things.
% Used in resfields*, resfreqs*, endorfrq*, curry, saffron, resfields_eig

%{
Inputs:
  Exp.PowderSimulation:   true/false (defined if coming from pepper/salt/saffron/curry)
  Exp.SampleFrame:        set by pepper/salt/saffron for powders, by user for cystals
  Exp.SampleRotation      rotation of the sample in the lab frame
  Exp.CrystalSymmetry:    space group, only used for crystal simulations
  Exp.MolFrame:           molecular frame orientation in crystal frame, only used for crystal sims
  Exp.R_sample            sample rotation matrix
  Opt.Sites:              list of sites to include, e.g. [1 3] in a 4-site group

Outputs:
  Orientations:    list of orientations Euler angles, one per row
                     [phi,theta,chi] for mol-to-lab frame transformation
  nOrientations:   total number of orientations
  nSites:          number of symmetry-related sites
  averageOverChi:  whether to compute average over third angle
%}

if ~isfield(Exp,'SampleFrame')
  error('Internal error: Exp.SampleFrame missing. Please report.');
end
% SampleFrame contains an Nx3 array of Euler angles (in radians)
% specifying the relative orientation between laboratory and sample frame.
if size(Exp.SampleFrame,2)~=3
  error('Exp.SampleFrame requires three columns (three Euler angles).')
end

rotateSample = isfield(Exp,'R_sample') && ~isempty(Exp.R_sample);

% Exp.PowderSimulation is set only if this function is called from an EasySpin
% function that does a powder simulation (pepper, salt, saffron, curry, etc).
doPowderSimulation = isfield(Exp,'PowderSimulation') && Exp.PowderSimulation;

if ~doPowderSimulation
  
  % Crystals
  %-----------------------------------------------------------------------------------
  % Get site-to-site transformations for given crystal space group
  % (matrices for active rotations; defined in crystal frame)
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
  
  % Crystal-to-molecular frame transformation, R_C2M
  if ~isempty(Exp.MolFrame)
    logmsg(1,'  molecular frame given in Exp.MolFrame');
    if ~numel(Exp.MolFrame)==3
      error('Exp.MolFrame has wrong size.');
    end
    R_C2M = erot(Exp.MolFrame);
  else
    R_C2M = eye(3);
  end
  R_M2C = R_C2M.';
  
  % Generate list of mol frame orientations, represented in crystal/sample frame
  xyzM_M = eye(3);           % xM, yM, zM in mol frame representation
  xyzM0_C = R_M2C*xyzM_M;   % transformation to crystal/sample frame
  for iSite = nSites:-1:1
    xyzM_C{iSite} = Rsite_C{iSite}*xyzM0_C;   % active rotation in crystal frame
  end

  % Lab-to-crystal/sample frame transformation, R_L2C
  % - R_L2C col 1,2,3: lab axis 1,2,3 represented in crystal frame
  % - R_L2C row 1,2,3: crystal axis 1,2,3 represented in lab frame
  if ~isnumeric(Exp.SampleFrame)
    error('Exp.SampleFrame must be an Nx3 array.');
  end
  if ~isempty(Exp.SampleFrame)
    nSamples = size(Exp.SampleFrame,1);
    logmsg(1,'  %d sample orientation(s) given in Exp.SampleFrame',nSamples);
    
    % Construct transformation matrices (lab frame to crystal frame)
    for s = nSamples:-1:1
      R_L2C{s} = erot(Exp.SampleFrame(s,:));
    end
    
  else
    R_L2C = {eye(3)};
  end
  nOrientations = numel(R_L2C);
  
  % Apply (active) sample rotation
  if rotateSample
    for iR = 1:numel(R_L2C)
      R_L2C{iR} = Exp.R_sample*R_L2C{iR};
    end
  end
  
  % Generate list of lab frame orientations, represented in the
  % various site molecular frames
  % (M = molecular frame of reference site)
  % (Mi = molecular frame of site # i)
  angles_M2L = zeros(nOrientations*nSites,3);
  idx = 1;
  for iOri = 1:nOrientations
    R_C2L = R_L2C{iOri}.';
    for iSite = 1:nSites
      xyzMi_L = R_C2L*xyzM_C{iSite};    % transformation to lab frame
      angles_M2L(idx,:) = eulang(xyzMi_L,true); % Euler angles for Mi -> L frame
      idx = idx + 1;
    end
  end
  nOrientations = nOrientations*nSites;
  
  averageOverChi = false;
  
else
  
  % Powder
  %-----------------------------------------------------------------------------------
  % Powder simulation: Orientations supplied by pepper etc in Exp.SampleFrame.
  
  % Exp.SampleFrame transforms from lab to sample frame, but Orientations = [phi,
  % theta, chi] from mol to lab frame. (MolFrame is ignored for a powder)
  angles_M2L = -fliplr(Exp.SampleFrame);
  nOrientations = size(angles_M2L,1);
  
  % For powder simulations, always average over the third angle.
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
