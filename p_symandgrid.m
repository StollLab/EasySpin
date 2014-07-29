%==========================================================================
% Symmetry determination and orientational grid.
%==========================================================================
logmsg(1,'-orientations------------------------------------------');

if (PowderSimulation)
  
  logmsg(1,'  powder sample (randomly oriented centers)');
  
  if isempty(Opt.Symmetry)
    
    msg = 'automatic determination of symmetry group and frame';
    [Opt.Symmetry,Opt.SymmFrame] = symm(Sys);
    if (numel(Exp.Temperature)>1) && strcmp(Opt.Symmetry,'Dinfh')
      logmsg(1,'  Hamiltonian symmetry is axial, non-equilibrium populations\n   -> reduction to rhombic');
      Opt.Symmetry = 'D2h';
    end
    if isfield(Opt,'ThetaRange')
      if ~isempty(Opt.ThetaRange)
        Opt.Symmetry = 'Ci';
        Opt.SymmFrame = eye(3);
      end
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
  [maxPhi,openPhi,nOctants] = symparam(Opt.Symmetry);
  openPhi = 0;
  InterpParams = [Opt.nKnots(1),openPhi,nOctants];
  
  
  % Display symmetry group and frame.
  if (TiltedFrame), msgg = 'tilted'; else msgg = 'non-tilted'; end
  logmsg(1,'  %s, %s frame',Opt.Symmetry,msgg);
  if (TiltedFrame)
    str = '  symmetry frame orientation in reference frame (Euler angles, deg):\n    [%s]';
    logmsg(1,str,sprintf('%0.3f ',eulang(Opt.SymmFrame,1)*180/pi));
  end
  
  % Get orientations for the knots, molecular frame.
  [Vecs,Weights] = sphgrid(Opt.Symmetry,Opt.nKnots(1),'cf');
  
  % Transform vector to reference frame representation and convert to polar angles.
  [phi,theta] = vec2ang(Opt.SymmFrame*Vecs);
  clear Vecs;
  Exp.CrystalOrientation = [phi;theta].';
  nOrientations = numel(phi);
  
  % Display information on orientational grid
  switch (nOctants)
    case -1, str = '  region: north pole (single point)';
    case 0,  str = '  region: a quarter of a meridian';
    otherwise, str = sprintf('  region: %d octant(s) of the upper hemisphere',nOctants);
  end
  logmsg(1,str);
  logmsg(1,'  %d orientations (%d knots)',nOrientations,Opt.nKnots(1));
  
else % no powder simulation
  
  nOrientations = size(Exp.CrystalOrientation,1);

  openPhi = 1;
  nOctants = -2;
  
  Weights = ones(1,nOrientations)*4*pi; % for consistency with powder sims (factor from integral over phi and theta)
  
  if (nOrientations==1)
    logmsg(1,'  crystal sample');
  else
    logmsg(1,'  crystal sample with %d user-specified orientations',nOrientations);
  end
  
end
%==========================================================================
