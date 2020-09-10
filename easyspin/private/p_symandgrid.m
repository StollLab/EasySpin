function [Exp,Opt] = p_symandgrid(Sys,Exp,Opt)
%==========================================================================
% Symmetry determination and orientational grid.
%==========================================================================
logmsg(1,'-orientations------------------------------------------');

if Exp.PowderSimulation
  
  logmsg(1,'  powder sample (randomly oriented centers)');
  
  if isempty(Opt.Symmetry)
    
    msg = 'automatic determination of symmetry group and frame';
    [Opt.Symmetry,Opt.SymmFrame] = symm(Sys);
    nonEquiPops = isfield(Sys,'Pop') && ~isempty(Sys.Pop);
    if nonEquiPops && strcmp(Opt.Symmetry,'Dinfh')
      logmsg(1,'  Hamiltonian symmetry is axial, non-equilibrium populations\n   -> reduction to rhombic');
      Opt.Symmetry = 'D2h';
    end
  
  else

    msg = 'user-specified symmetry group';
    if isempty(Opt.SymmFrame)
      Opt.SymmFrame = eye(3);
      msg = [msg ', fixed symmetry frame'];
    else
      msg = [msg ' and symmetry frame'];
    end
  end
  logmsg(1,'  %s',msg);

  TiltedFrame = (sum(diag(Opt.SymmFrame))~=3);
  
  
  % Convert symmetry to octant number.
  [maxPhi,Opt.openPhi,Opt.nOctants] = symparam(Opt.Symmetry);
  Opt.openPhi = 0;
  Opt.InterpParams = [Opt.nKnots(1),Opt.openPhi,Opt.nOctants];
  
  
  % Display symmetry group and frame.
  if (TiltedFrame), msgg = 'tilted'; else, msgg = 'non-tilted'; end
  logmsg(1,'  %s, %s frame',Opt.Symmetry,msgg);
  if (TiltedFrame)
    str = '  symmetry frame orientation in reference frame (Euler angles, deg):\n    [%s]';
    logmsg(1,str,sprintf('%0.3f ',eulang(Opt.SymmFrame,1)*180/pi));
  end
  
  % Get orientations for the knots, molecular frame.
  [Vecs,Exp.OriWeights] = sphgrid(Opt.Symmetry,Opt.nKnots(1),'cf');
  
  % Transform vector to reference frame representation and convert to polar angles.
  [Exp.phi,Exp.theta] = vec2ang(Opt.SymmFrame*Vecs);
  clear Vecs
  Exp.CrystalOrientation = [Exp.phi;Exp.theta].';
  nOrientations = numel(Exp.phi);
  
  % Display information on orientational grid
  switch Opt.nOctants
    case -1, str = '  region: north pole (single point)';
    case 0,  str = '  region: a quarter of a meridian';
    otherwise, str = sprintf('  region: %d octant(s) of the upper hemisphere',Opt.nOctants);
  end
  logmsg(1,str);
  if Opt.nOctants==-1
    logmsg(1,'  1 orientation (1 knot)');
  else
    logmsg(1,'  %d orientations (%d knots)',nOrientations,Opt.nKnots(1));
  end
  
else % no powder simulation

  % Transpose Exp.CrystalOrientation if necessary to have one row per orientation
  transpose = false;
  nC1 = size(Exp.CrystalOrientation,1);
  nC2 = size(Exp.CrystalOrientation,2);
  if nC2==2 || nC2==3
    % all fine
  else
    if nC1==2 || nC1==3
      transpose = true;
    else
      error('Exp.CrystalOrientation must be a Nx3 or Nx2 array, yours is %dx%d.',...
        nC1,nC2);
    end
  end
  if transpose, Exp.CrystalOrientation = Exp.CrystalOrientation.'; end
  nOrientations = size(Exp.CrystalOrientation,1);

  Opt.openPhi = 1;
  Opt.nOctants = -2;
  
  Exp.OriWeights = ones(1,nOrientations)*4*pi; % for consistency with powder sims (factor from integral over phi and theta)
  
  if nOrientations==1
    logmsg(1,'  crystal sample');
  else
    logmsg(1,'  crystal sample with %d user-specified orientations',nOrientations);
  end
  
end
