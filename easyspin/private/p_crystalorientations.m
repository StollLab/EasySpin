function [Orientations,nOrientations,nSites,averageOverChi] = p_crystalorientations(Exp,Opt)

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
  
  % Crystal simulation
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
  
  % Crystal-to-molecular frame transformation, R_CM
  % - R_CM col 1,2,3: crystal axis xC,yC,zC represented in molecular frame
  % - R_CM row 1,2,3: molecular axis xM,yM,zM represented in crystal frame
  if ~isempty(Exp.MolFrame)
    logmsg(1,'  molecular frame given in Exp.MolFrame');
    if isequal(size(Exp.MolFrame),[3 3])
      R_CM = Exp.MolFrame;
    elseif numel(Exp.MolFrame)==3
      R_CM = erot(Exp.MolFrame);
    else
      error('Exp.MolFrame has wrong size.');
    end
  else
    R_CM = eye(3);
  end
  
  % Crystal-to-lab frame transformation, R_CL
  % - R_CL col 1,2,3: crystal axis 1,2,3 represented in lab frame
  % - R_CL row 1,2,3: lab axis 1,2,3 represented in crystal frame
  if ~isnumeric(Exp.SampleFrame)
    error('Exp.SampleFrame must be an Nx3 array.');
  end
  if ~isempty(Exp.SampleFrame)
    nSamples = size(Exp.SampleFrame,1);
    logmsg(1,'  %d sample orientation(s) given in Exp.SampleFrame',nSamples);
    
    % Construct transformation matrices (lab frame to crystal frame)
    for s = nSamples:-1:1
      R_LC{s} = erot(Exp.SampleFrame(s,:));
    end
    
  else
    R_LC = {eye(3)};
  end
  
  nOrientations = numel(R_LC);
  
  % Apply (active) sample rotation
  if rotateSample
    for iR = 1:numel(R_LC)
      R_LC{iR} = Exp.R_sample*R_LC{iR};
    end
  end
  
  % Generate list of lab frame orientations, represented in the
  % various site molecular frames
  % (M = molecular frame of reference site)
  % (Mi = molecular frame of site # i)
  xyzM_M = eye(3);  % xM, yM, zM in mol frame representation
  xyzM_C = R_CM.'*xyzM_M;
  Orientations = zeros(nOrientations*nSites,3);
  idx = 1;
  for iOri = 1:nOrientations
    for iSite = 1:nSites
      xyzMi_C = Rsite_C{iSite}*xyzM_C;
      xyzMi_L = R_LC{iOri}.'*xyzMi_C;
      xyzL_Mi = xyzMi_L.';
      Orientations(idx,:) = eulang(xyzL_Mi.',true);
      idx = idx + 1;
    end
  end
  nOrientations = nOrientations*nSites;
  
  averageOverChi = false;
  
else
  
  % Powder simulation
  %-----------------------------------------------------------------------------------
  % Powder simulation: Orientations supplied by pepper etc in Exp.SampleFrame.
  
  % Exp.SampleFrame transforms from lab to sample frame, but Orientations = [phi,
  % theta, chi] from mol to lab frame. (MolFrame is ignored for a powder)
  Orientations = -fliplr(Exp.SampleFrame);
  nOrientations = size(Orientations,1);
  
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
