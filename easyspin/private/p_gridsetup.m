% Symmetry determination and orientational grid
%
% requires
%   Sys (for hamsymm() call)
%   Sys.initState (optional)
%   Exp.SampleFrame (for crystals)
%   Opt.GridIntegration (required)
%   Opt.GridSymmetry (optional)
%   Opt.GridFrame (optional)
%   Opt.GridSize
%
% sets
%   for all samples
%     Exp.OriWeights
%     Opt.nOctants
%   for disordered samples:
%     Exp.MolFrame
%     Exp.phi
%     Exp.theta
%     Exp.tri
%     Opt.GridSymmetry
%     Opt.GridFrame
%     Opt.GridParams

function [Exp,Opt,nOrientations] = p_gridsetup(Sys,Exp,Opt)

logmsg(1,'-orientations------------------------------------------');

if Opt.GridIntegration
  % disordered sample (powder or frozen solution) or partially ordered sample
  
  logmsg(1,'  disordered sample (randomly oriented centers)');
  
  if isempty(Opt.GridSymmetry)
    
    % Determine grid symmetry from Hamiltonian symmetry
    msg = 'automatic determination of grid symmetry and orientation';
    [Opt.GridSymmetry,Opt.GridFrame] = hamsymm(Sys);
    nonEquiPops = isfield(Sys,'initState') && ~isempty(Sys.initState);
    if nonEquiPops && strcmp(Opt.GridSymmetry,'Dinfh')
      logmsg(1,'  Hamiltonian symmetry is axial (Dinfh), non-equilibrium state given in Sys.initState\n   -> reducing symmetry to rhombic (D2h).');
      Opt.GridSymmetry = 'D2h';
    end
  
  else
    
    msg = 'user-specified grid symmetry';
    if isempty(Opt.GridFrame)
      Opt.GridFrame = eye(3);
    else
      msg = [msg ' and grid frame orientation'];
    end
  end
  logmsg(1,'  %s',msg);
  
  tiltedFrame = sum(diag(Opt.GridFrame))~=3; % test for eye(3)
  
  % Display grid symmetry group and grid frame
  if tiltedFrame
    msg = 'tilted';
  else
    msg = 'non-tilted';
  end
  logmsg(1,'  grid: symmetry %s, frame %s',Opt.GridSymmetry,msg);
  if tiltedFrame
    str = '  grid frame orientation in molecular frame (Euler angles, deg):\n    [%s]';
    logmsg(1,str,sprintf('%0.3f ',eulang(Opt.GridFrame,1)*180/pi));
  end
  
  % Get spherical grid
  [grid,tri] = sphgrid(Opt.GridSymmetry,Opt.GridSize(1));
  Vecs = grid.vecs;
  Opt.nOctants = grid.nOctants;
  Exp.OriWeights = grid.weights;
  Exp.tri = tri;
  
  % Transform grid vectors to molecular frame representation, convert to
  % polar angles, and store in Exp.MolFrame
  [Exp.phi,Exp.theta] = vec2ang(Opt.GridFrame*Vecs);
  chi = zeros(size(Exp.theta));
  Exp.MolFrame = [-chi; -Exp.theta; -Exp.phi].';
  nOrientations = numel(Exp.phi);
  
  Opt.GridParams = [nOrientations,Opt.GridSize(1),grid.closedPhi,Opt.nOctants,grid.maxPhi];

  % Display information about orientational grid
  switch grid.nOctants
    case -1, str = 'north pole (single point)';
    case 0, str = 'a quarter of a meridian';
    case 4, str = 'upper hemisphere';
    case 2, str = 'half of upper hemisphere';
    case 8, str = 'full sphere';
    case 1, str = sprintf('%d octant of the upper hemisphere, maxPhi = %g*pi',Opt.nOctants,grid.maxPhi/pi);
    otherwise
      error('Unsupported number %d of grid octants.',Opt.nOctants);
  end
  logmsg(1,['  integration region: ' str]);
  if nOrientations==1
    logmsg(1,'  1 orientation (GridSize = %d)',Opt.GridSize(1));
  else
    logmsg(1,'  %d orientations (GridSize = %d)',nOrientations,Opt.GridSize(1));
  end
  
else  % crystal sample 
  
  % Check Exp.SampleFrame
  [nC1,nC2] = size(Exp.SampleFrame);
  if nC2~=2 && nC2~=3
    error('Exp.SampleFrame must be a Nx3 or Nx2 array, yours is %dx%d.',...
        nC1,nC2);
  end
  nOrientations = size(Exp.SampleFrame,1);
  
  Opt.nOctants = -2;
  
  Exp.OriWeights = ones(1,nOrientations)*4*pi;
  % for consistency with powder sims (factor from integral over phi and theta)
  
  if nOrientations==1
    logmsg(1,'  crystal sample, single orientation');
  else
    logmsg(1,'  crystal sample with %d user-specified orientations',nOrientations);
  end
  
end

