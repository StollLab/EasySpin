% Script that takes care of crystal frame setups for single-crystal
% simulation and a few related things.
% Used in resfields*, resfreqs*, endorfrq*, curry

%{
inputs:
  Exp.PowderSimulation:   true/false if coming from pepper/salt/saffron/curry, not defined otherwise
  Exp.CrystalOrientation: set by pepper/salt/saffron for powders, by user for cystals
  Exp.CrystalSymmetry
  Exp.MolFrame
  Opt.Sites

outputs:
  Orientations:    list of orientations Euler angles, one per row
  nSites:          number of symmetry-related sites
  AverageOverChi:  whether to compute average over third angle
%}

if isfield(Exp,'Orientations')
  error('Exp.Orientations is obsolete. Use Exp.CrystalOrientation instead.');
end

% Exp.PowderSimulation is set only if this is an internal call via
% pepper, salt or saffron.
if ~isfield(Exp,'PowderSimulation')
  PowderSimulation = false;
else
  PowderSimulation = Exp.PowderSimulation;
end

if (~PowderSimulation)
    
  % Get site-to-site transformations for given crystal space group
  % (matrices for active rotations; defined in crystal frame)
  if isempty(Exp.CrystalSymmetry)
    Exp.CrystalSymmetry = 'P1';
    logmsg(1,sprintf('  no crystal symmetry given: assuming space group P1'));
  end
  Rsite_C = sitetransforms(Exp.CrystalSymmetry);
  nSites  = numel(Rsite_C);
  logmsg(1,sprintf('  crystal symmetry: %d magnetically distinct sites in unit cell',nSites));
  
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
    logmsg(1,'  molecular frame given');
    if all(size(Exp.MolFrame)==[3 3])
      R_CM = Exp.MolFrame;
    elseif numel(Exp.MolFrame)==3
      R_CM = erot(Exp.MolFrame);
    else
      error('Exp.MolFrame has wrong size.');
    end
  else
    R_CM = 1;
  end
  
  % Crystal-to-lab frame transformation, R_CL
  % - R_CL col 1,2,3: crystal axis 1,2,3 represented in lab frame
  % - R_CL row 1,2,3: lab axis 1,2,3 represented in crystal frame
  if ~isempty(Exp.CrystalOrientation)
    logmsg(1,'  crystal orientation(s) given');

    % Transpose Exp.CrystalOrientation if necessary to have
    % one row per orientation
    transpose = false;
    nC1 = size(Exp.CrystalOrientation,1);
    nC2 = size(Exp.CrystalOrientation,2);
    if (nC2==2) || (nC2==3)
      % all fine
    else
      if (nC1==2) || (nC1==3)
        transpose = true;
      else
        error('Exp.CrystalOrientation must be a Nx3 or Nx2 array, yours is %dx%d.',...
          nC1,nC2);
      end
    end
    if transpose, Exp.CrystalOrientation = Exp.CrystalOrientation.'; end
        
    COri = Exp.CrystalOrientation;
    if all(size(COri)==[3 3])
      R_CL = COri;
    elseif numel(COri)==2
      R_CL = erot([COri(1) COri(2) 0]);
    elseif numel(COri)==3
      R_CL = erot(COri);
    elseif any(size(COri)==3)
      if size(COri,1)==3, COri = COri.'; end
      for iOri = 1:size(COri,1)
        R_CL{iOri} = erot(COri(iOri,:));
      end
    elseif any(size(COri)==2)
      if size(COri,1)==2, COri = COri.'; end
      for iOri = 1:size(COri,1)
        R_CL{iOri} = erot([COri(iOri,:) 0]);
      end
    else
      error('Exp.CrystalOrientation has wrong size.');
    end
    Exp.CrystalOrientation = COri;
  else
    R_CL = 1;
  end
  if ~iscell(R_CL), R_CL = {R_CL}; end
  nOrientations = numel(R_CL);
  
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
  
  AverageOverChi = false;

else
  
  if isfield(Exp,'CrystalOrientation')
    % Powder simulation: Orientations supplied by pepper in
    % Exp.CrystalOrientations, Exp.CrystalSymmetry and Exp.MolFrame
    % are disregarded, since they have no effect on the powder spectrum.
    nSites = 1;
    Orientations = Exp.CrystalOrientation;
    
    if ~isempty(Exp.CrystalSymmetry)
      logmsg(1,'  powder simulation; disregarding Exp.CrystalSymmetry');
    end
    if ~isempty(Exp.MolFrame)
      logmsg(1,'  powder simulation; disregarding Exp.MolFrame');
    end
    
    % Contains a nx2 or nx3 array of angles (in radian units). These specify
    % the relative orientation between molecular and laboratory frame
    % [phi theta chi]. If the third angle is missing, it is set to zero.
    [nOrientations,nAngles] = size(Orientations);
    switch nAngles
      case 2
        Orientations(end,3) = 0; % Entire chi row is set to 0.
      case 3
        % ok
      otherwise
        error('Orientations array has %d columns instead of 2 or 3.',nAngles);
    end
  else
    error('Internal error. Please report.');
  end
  
  AverageOverChi = true;
  
end
