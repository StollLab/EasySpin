% orisel  Orientation selection 
%
%   orisel(Sys,Exp)
%   orisel(Sys,Exp,Opt)
%   Weights = ...
%   [Weights,Trans] = ...
%
%   Calculates weights for orientation selectivity.
%
%   Input:
%   - Sys    spin system
%   - Exp    experimental parameters
%      Field         magnetic field [mT]
%      mwFreq        spectrometer frequency [GHz]
%      ExciteWidth   FWHM excitation bandwidth [MHz]
%      Orientations  [phi;theta]
%   - Opt    options structure
%      GridSize      grid dsize
%      GridSymmetry   grid symmetry: C1, Ci, D2h, or Dinfh
%
%   Output:
%   - Weights   vector of selectivity weights, one per orientation
%   - Trans     list of transitions included
%
%   If the function is called with no output arguments,
%   the selectivity weights are plotted.

% Syntax with 2 output arguments is undocumented since 2.2.0

function varargout = orisel(Sys,Params,Options)

if nargin==0, help(mfilename); return; end

error(chkmlver);

if nargin<3, Options = struct; end

% Process spin system.
%----------------------------------------------------------------------
[Sys,err] = validatespinsys(Sys);
error(err);
if ~isfield(Sys,'HStrain'), Sys.HStrain = [0 0 0]; end
if numel(Sys.HStrain)==1, Sys.HStrain = [1 1 1]*Sys.HStrain; end

% Process Parameters.
%----------------------------------------------------------------------
DefaultParams.mwFreq = NaN;
DefaultParams.Field = NaN;
DefaultParams.ExciteWidth = inf;
DefaultParams.Orientations = [];

if ~isfield(Params,'Field'), error('Params.Field is missing!'); end
if ~isfield(Params,'mwFreq'), error('Params.mwFreq is missing!'); end
if ~isfield(Params,'ExciteWidth'), error('Params.ExciteWidth is missing!'); end

Params = adddefaults(Params,DefaultParams);

if any(Sys.HStrain+Params.ExciteWidth<=0)
  error('Sys.HStrain or Exp.ExciteWidth must be positive!');
end

mwFreq = Params.mwFreq*1e3; % GHz -> MHz

% Process Options
%-----------------------------------------------------------------------
DefaultOptions.GridSize = 46; % same as in salt and nucfrq2d
DefaultOptions.GridSymmetry = '';
DefaultOptions.Display = (nargout==0);

DefaultOptions.AveragedIntensity = true;

Options = adddefaults(Options,DefaultOptions);

if ~isempty(Params.Orientations)
  Orientations = Params.Orientations;
  phi = Orientations(1,:);
  theta = Orientations(2,:);
else
  if isempty(Options.GridSymmetry)
    [Options.GridSymmetry,GridFrame] = symm(Sys);
  else
    GridFrame = eye(3);
  end
  SymmGroup = Options.GridSymmetry;
  if Options.Display
    % no open phi intervals, since triangulation only works for closed ones
    gridopt = 'c'; 
  else
    gridopt = '';
  end
  [grid,tri] = sphgrid(SymmGroup,Options.GridSize,gridopt);
  Vectors = grid.vecs;
  [phi,theta] = vec2ang(GridFrame*Vectors);
  Orientations = [phi; theta];
end

[nAngles,nOrientations] = size(Orientations);
switch nAngles
  case 2
    Orientations(3,end) = 0;
  case 3
    % don't do anything
  otherwise
    error('Orientations array has wrong size, should be 2xn or 3xn.');
end


%-----------------------------------------------------------------------
% Prepare Hamiltonian and get state space dimension
[F,GxM,GyM,GzM] = sham(Sys);
N = length(F);

% Transition map and level indices
uv = find(triu(ones(N),1));
[u,v] = ind2sub([N N],uv); % u < v
Transitions = [u,v];
nTransitions = numel(u);

% Pre-calculate
fac = 1/sqrt(2*log(2));

if isfinite(Params.ExciteWidth)

  Weights = zeros(nOrientations,nTransitions);
  % Pre-calculate width info
  exc2 = Params.ExciteWidth^2;
  HStrain2 = Sys.HStrain.^2;
  one = ones(1,nOrientations);
  vec = ang2vec(Orientations(1,:),Orientations(2,:)).';
  GammaWidth = fac*sqrt(sum(vec.^2.*HStrain2(one,:),2) + exc2(one,1));
  
  % Calculate weights for all orientations and all transitions
  for iOri = 1:nOrientations
    [xLab,yLab,zLab] = erot(phi(iOri),theta(iOri),0,'rows');
    GxL = xLab(1)*GxM + xLab(2)*GyM + xLab(3)*GzM;
    GyL = yLab(1)*GxM + yLab(2)*GyM + yLab(3)*GzM;
    GzL = zLab(1)*GxM + zLab(2)*GyM + zLab(3)*GzM;
    % Eigenvalues
    [V,E] = eig(F + Params.Field*GzL);
    E = diag(E);
    [E,idx] = sort(E); % because of a bug in eig() in Matlab 7.0.0 (fixed in 7.0.1)
    V = V(idx,:);
    
    if (Options.AveragedIntensity)
      TransitionRate = (abs(V'*GxL*V).^2 + abs(V'*GyL*V).^2)/2;
    else
      TransitionRate = abs(V'*kGxL*V).^2;
    end

    % Orientation selectivity weights
    xi = (abs(E(v)-E(u))-mwFreq)/GammaWidth(iOri);
    Weights(iOri,:) = (TransitionRate(uv) .* exp(-2*xi.^2)).';
  end

else
  
  Weights = ones(nOrientations,nTransitions);

end

if Options.Display
  if isempty(Params.Orientations)
    if strcmp(SymmGroup,'Dinfh')
      plot(theta*180/pi,sum(Weights,2));
      xlabel('theta [deg]');
      ylabel('weights [a.u.]');
    else
      trisurf(tri.idx,vec(:,1),vec(:,2),vec(:,3),sum(Weights,2));
      view([130 40]);
      c = caxis; c(1) = 0; caxis(c);
      shading flat
      axis equal
      colorbar;
      xlabel('x_{M}');
      ylabel('y_{M}');
      zlabel('z_{M}');
    end
  else
    disp('Cannot plot orientation selection for a user-supplied set of orientations.');
  end
end

switch nargout
  case 1, varargout = {sum(Weights,2)};
  case 2, varargout = {Weights,Transitions};
end

return
