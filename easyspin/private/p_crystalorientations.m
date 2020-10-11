function [Orientations,nOrientations,nSites,averageOverChi] = p_crystalorientations(Exp,Opt)

% Function that takes care of crystal frame setups for single-crystal simulation
%    and a few related things.
% Used in resfields*, resfreqs*, endorfrq*, curry, saffron, eigfields

%{
Inputs:
  Exp.PowderSimulation:   true/false (defined if coming from pepper/salt/saffron/curry)
  Exp.CrystalOrientation: set by pepper/salt/saffron for powders, by user for cystals
  Exp.SampleRotation      rotation of the sample in the lab frame
  Exp.CrystalSymmetry:    space group, only used for crystal simulations
  Exp.MolFrame:           molecular frame orientation in crystal frame, only used for crystal sims
  Opt.Sites:              list of sites to include, e.g. [1 3] in a 4-site group

Outputs:
  Orientations:    list of orientations Euler angles, one per row
  nOrientations:   total number of orientations
  nSites:          number of symmetry-related sites
  AverageOverChi:  whether to compute average over third angle
%}

% Exp.PowderSimulation is set only if this function is called from an EasySpin
% function that does a powder simulation (pepper, salt, saffron, curry, etc).
PowderSimulation = isfield(Exp,'PowderSimulation') && Exp.PowderSimulation;

% If sample rotation is requested, prepare rotation matrices
rotateSample = isfield(Exp,'SampleRotation') && ~isempty(Exp.SampleRotation);
if rotateSample
  if isnumeric(Exp.SampleRotation)
    rho = Exp.SampleRotation;
    nRot = [1;0;0]; % lab x axis (xL_L) as default
  elseif iscell(Exp.SampleRotation) && numel(Exp.SampleRotation)==2
    rho = Exp.SampleRotation{1};
    nRot = Exp.SampleRotation{2};
  else
    error('Exp.SampleRotation must be of the form {rho,nL}, with the rotation angle rho and the lab-frame rotation axis nL.');
  end
  if numel(rho)~=1
    error('Exp.SampleRotation: the first element must the rotation angle rho.');
  end
  R_sample = rotaxi2mat(nRot,rho).'; % transpose because it's an active rotation
end

if ~PowderSimulation
  
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
    if all(size(Exp.MolFrame)==[3 3])
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
  if ~isnumeric(Exp.CrystalOrientation)
    error('Exp.CrystalOrientation must be a Nx2 or Nx3 array.');
  end
  if ~isempty(Exp.CrystalOrientation)
    logmsg(1,'  crystal orientation(s) given in Exp.CrystalOrientation');
    
    % Check Exp.CrystalOrientation, supplement third angle if necessary.
    Ori_ = Exp.CrystalOrientation;
    nC1 = size(Ori_,1);
    nC2 = size(Ori_,2);
    if nC2==3
      % Nx3: ok, do nothing
    elseif nC2==2
      % Nx2: supplement third angle chi
      Ori_(:,3) = 0;
    else
      error('Exp.CrystalOrientation must be a Nx2 or Nx3 array, yours is %dx%d.',nC1,nC2);
    end
    Exp.CrystalOrientation = Ori_;
    
    % Construct rotation matrices (crystal frame to lab frame)
    for iOri = size(Ori_,1):-1:1
      R_CL{iOri} = erot(Ori_(iOri,:));
    end
    
  else
    R_CL = {eye(3)};
  end
  
  nOrientations = numel(R_CL);
  
  % Apply sample rotation
  if rotateSample
    for iR = 1:numel(R_CL)
      R_CL{iR} = R_sample*R_CL{iR};
    end
  end
  
  % Generate list of lab frame orientations, represented in the
  % various site molecular frames
  % (M = molecular frame of reference site)
  % (Mi = molecular frameof site # i)
  xyzM_M = eye(3);
  xyzM_C = R_CM.'*xyzM_M;    
  Orientations = zeros(nOrientations*nSites,3);
  idx = 1;
  for iOri = 1:nOrientations
    for iSite  = 1:nSites
      xyzMi_C = Rsite_C{iSite}*xyzM_C;
      xyzMi_L = R_CL{iOri}*xyzMi_C;
      xyzL_Mi = xyzMi_L.';
      Orientations(idx,:) = eulang(xyzL_Mi.',1);
      idx = idx + 1;
    end
  end
  nOrientations = nOrientations*nSites;
  
  averageOverChi = false;
  
else
  
  % Powder simulation
  %-----------------------------------------------------------------------------------
  % Powder simulation: Orientations supplied by pepper etc in Exp.CrystalOrientation.
  if ~isfield(Exp,'CrystalOrientation')
    error('Internal error: No orientations for powder simulation. Please report.');
  end
  
  % Orientations contains a nx2 or nx3 array of Euler angles (in radian units).
  % These specify the relative orientation between molecular and laboratory frame
  % [phi theta chi]. If the third angle, chi, is missing, it is set to zero.
  Orientations = Exp.CrystalOrientation;
  [nOrientations,nAngles] = size(Orientations);
  switch nAngles
    case 2
      Orientations(end,3) = 0; % Entire chi row is set to 0.
    case 3
      % ok
    otherwise
      error('Orientations array has %d columns instead of 2 or 3.',nAngles);
  end
  
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
