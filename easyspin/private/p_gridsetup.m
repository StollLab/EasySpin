% Symmetry determination and orientational grid
%
% requires
%   Sys (for hamsymm() call)
%   Sys.initState (optional)
%   Exp.SampleFrame (for crystals)
%   Opt.disorderedSample (required)
%   Opt.partiallyOrderedSample (required)
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

logmsg(1,'-grid/orientations-------------------------------------');

useGrid = Opt.disorderedSample || Opt.partiallyOrderedSample;
if useGrid
  % disordered sample (powder or frozen solution)
  
  logmsg(1,'  disordered sample (randomly oriented centers)');
  
  % Set grid symmetry group
  if isempty(Opt.GridSymmetry)    
    if Opt.disorderedSample
      % Determine grid symmetry from Hamiltonian symmetry
      [Opt.GridSymmetry,Opt.GridFrame] = hamsymm(Sys);
      nonEquiPops = isfield(Sys,'initState') && ~isempty(Sys.initState);
      if nonEquiPops && strcmp(Opt.GridSymmetry,'Dinfh')
        logmsg(1,'  Hamiltonian symmetry is axial (Dinfh), non-equilibrium state given in Sys.initState\n   -> reducing symmetry to rhombic (D2h).');
        Opt.GridSymmetry = 'D2h';
      end
    else
      % Use Ci grid symmetry for partially ordered samples (since symmetry
      % of ordering function is unknown)
      if isempty(Opt.GridSymmetry)
        Opt.GridSymmetry = 'Ci';
      end
    end
    msg = '  grid: automatically determined grid symmetry';
  else
    msg = '  grid: user-specified grid symmetry';
  end
  logmsg(1,msg);

  if Opt.partiallyOrderedSample && ~strcmp(Opt.GridSymmetry,'Ci')
    error('For partially ordered samples, Ci grid symmetry is required.');
  end

  % Set grid frame orientation
  if isempty(Opt.GridFrame)
    Opt.GridFrame = eye(3);
    msg = '  grid: default grid frame orientation';
  else
    msg = '  grid: user-specified grid frame orientation';
  end
  logmsg(1,msg);
  
  % Display grid symmetry group and grid frame orientation
  tiltedFrame = sum(diag(Opt.GridFrame))~=3; % test for eye(3)
  if tiltedFrame
    msg = 'tilted';
  else
    msg = 'non-tilted';
  end
  logmsg(1,'  grid: symmetry %s, %s frame',Opt.GridSymmetry,msg);
  if tiltedFrame
    msg = '  grid: frame orientation in molecular frame (Euler angles, deg):\n    [%s]';
    logmsg(1,msg,sprintf('%0.3f ',eulang(Opt.GridFrame,1)*180/pi));
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
  
else  % Crystal sample 
  
  % Check Exp.SampleFrame
  [nC1,nC2] = size(Exp.SampleFrame);
  if nC2~=3
    error('Exp.SampleFrame must be a Nx3 array, yours is %dx%d.',nC1,nC2);
  end
  nOrientations = size(Exp.SampleFrame,1);
  if isfield(Exp,'SampleRotation') && ~isempty(Exp.SampleRotation)
    nOrientations = nOrientations*numel(Exp.SampleRotation{2});
  end
  
  Opt.nOctants = -2;
  
  Exp.OriWeights = ones(1,nOrientations)*4*pi;
  % for consistency with powder sims (factor from integral over phi and theta)
  
  if nOrientations==1
    logmsg(1,'  crystal sample, single orientation');
  else
    logmsg(1,'  crystal sample, %d orientations',nOrientations);
  end
  
end
