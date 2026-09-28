% resfreqs_perturb Compute resonance frequencies for frequency-swept EPR
%
%   ... = resfreqs_perturb(Sys,Exp)
%   ... = resfreqs_perturb(Sys,Exp,Opt)
%   [Pos,Int] = resfreqs_perturb(...)
%   [Pos,Int,Wid] = resfreqs_perturb(...)
%   [Pos,Int,Wid,Trans] = resfreqs_perturb(...)
%
%   Computes frequency-domain EPR line positions, intensities and widths using
%   perturbation theory. Only systems with one electron spin are supported.
%
%   Input:
%    Sys: spin system structure
%    Exp: experimental parameters
%      Field               static field, in mT
%      Temperature         temperature, in K; if omitted: high-temperature limit
%      SampleFrame         Nx3 array of Euler angles (in radians) for sample/crystal orientations
%      CrystalSymmetry     crystal symmetry (space group etc.)
%      MolFrame            Euler angles (in radians) for molecular frame orientation
%      SampleRotation      sample rotation, {nRot_L rho}
%      mwMode              excitation mode: 'perpendicular', 'parallel', {k pol}
%                            pol: polarization angle, 'circular+', 'circular-', 'unpolarized'
%      lightBeam           photoexcitation: '', 'perpendicular', 'parallel', 'unpolarized', {k alpha}
%      lightScatter        isotropic fraction of photoexcitation, between 0 and 1
%    Opt: additional computational options
%      Verbosity           level of detail of printing; 0, 1, 2
%      PerturbOrder        perturbation order; 1 or 2
%      Sites               list of crystal sites to include (default []: all)
%
%   Output:
%    Pos     line positions (in MHz)
%    Int     line intensities
%    Wid     line widths, full width at half maximum (FWHM), in MHz
%    Trans   list of transitions (level indices; nuclear sublevel indices are approximate)

function varargout = resfreqs_perturb(Sys,Exp,Opt)

% Compute resonance frequencies based on formulas from Iwasaki, J.Magn.Reson. 16, 417-423 (1974)

% Assert correct Matlab version
warning(chkmlver);

% Check number of input arguments.
switch nargin
  case 0, help(mfilename); return;
  case 2, Opt = struct;
  case 3
  otherwise
    error('Use two or three inputs: resfreqs_perturb(Sys,Exp) or resfreqs_perturb(Sys,Exp,Opt)!');
end

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
logmsg(Opt.Verbosity);

% Spin system
%------------------------------------------------------
[Sys,err] = validatespinsys(Sys);
error(err);
S = Sys.S;
highSpin = any(S>1/2);

err = '';
if Sys.nElectrons~=1
  err = sprintf('Perturbation theory available only for systems with 1 electron. Yours has %d.',Sys.nElectrons);
end
if any(Sys.L(:))
  err = sprintf('Perturbation theory not available for electron spin coupled to orbital angular momentum!');
end
if any(Sys.DStrain(:)) && any(Sys.DFrame(:))
  err = 'D strain cannot be used with tilted D tensors.';
end
if any(strncmp(fieldnames(Sys),'Ham',3))
  err = 'Perturbation theory not available for higher order terms';
end
if isfield(Sys,'nn') && any(Sys.nn(:)~=0)
  err = 'Perturbation theory not available for nuclear-nuclear couplings (Sys.nn).';
end
if ~isempty(Sys.initState)
  err = 'Sys.initState is not supported by resfreqs_perturb.';
end
error(err);

if ~isfield(Sys,'gAStrainCorr')
  Sys.gAStrainCorr = +1;
else
  Sys.gAStrainCorr = sign(Sys.gAStrainCorr(1));
end

if Sys.fullg
  g = Sys.g;
else
  R_g2M = erot(Sys.gFrame).'; % g frame -> molecular frame
  g = R_g2M*diag(Sys.g)*R_g2M.';
end

if highSpin
  if Sys.fullD
    D = Sys.D;
  else
    R_D2M = erot(Sys.DFrame).'; % D frame -> molecular frame
    D = R_D2M*diag(Sys.D)*R_D2M.';
  end
  % make D traceless (required for Iwasaki expressions)
  D = D - eye(3)*trace(D)/3;
  trDD = trace(D^2);
end

nTransitions = 2*S; % number of allowed transitions for one electron spin

I = Sys.I;
nNuclei = Sys.nNuclei;
for iNuc = nNuclei:-1:1
  if Sys.fullA
    A_ = Sys.A((iNuc-1)*3+(1:3),:);
  else
    R_A2M = erot(Sys.AFrame(iNuc,:)).'; % A frame -> molecular frame
    A_ = R_A2M*diag(Sys.A(iNuc,:))*R_A2M.';
  end
  A{iNuc} = A_;
  detA(iNuc) = det(A_);
  if detA(iNuc)==0
    error('All hyperfine principal values must be non-zero.');
  end
  invA{iNuc} = inv(A_);
  trAA(iNuc) = trace(A_.'*A_);
  mI{iNuc} = -I(iNuc):I(iNuc);
end

% Table of mI values for all nuclear sublevels, one column per nucleus
if nNuclei>0
  mIc = allcombinations(mI{:});
  II1 = I.*(I+1);
else
  mIc = zeros(1,0);
end
nNucSublevels = size(mIc,1);


% Experiment
%------------------------------------------------------
DefaultExp.Field = NaN;
DefaultExp.Temperature = NaN;
DefaultExp.mwMode = 'perpendicular';

DefaultExp.SampleFrame = [0 0 0];
DefaultExp.CrystalSymmetry = '';
DefaultExp.MolFrame = [0 0 0];
DefaultExp.SampleRotation = [];

DefaultExp.lightBeam = '';
DefaultExp.lightScatter = 0;

Exp = adddefaults(Exp,DefaultExp);

% Check for obsolete fields
if isfield(Exp,'CrystalOrientation')
  error('Exp.CrystalOrientation is no longer supported, use Exp.SampleFrame/Exp.MolFrame instead.');
end
if isfield(Exp,'Mode')
  error('Exp.Mode is no longer supported. Use Exp.mwMode instead.');
end

err = '';
if ~isnumeric(Exp.Field) || numel(Exp.Field)~=1 || isnan(Exp.Field)
  err = 'Exp.Field (in mT) is missing or not a single number.';
elseif Exp.Field<0
  err = 'Exp.Field cannot be negative. Negative fields are only supported for field sweeps.';
end

error(err);

[xi1,xik,nB1,nk,nB0_L,mwmode] = p_excitationgeometry(Exp.mwMode);

useTemperature = p_temperature(Exp);

if ~isfield(Opt,'Sites')
  Opt.Sites = [];
end
if ~isfield(Opt,'separateSites'), Opt.separateSites = true; end % internal; false keeps sites as columns


% Photoselection
usePhotoSelection = ~isempty(Exp.lightBeam) && Exp.lightScatter<1;

if usePhotoSelection
  if ~isfield(Sys,'tdm') || isempty(Sys.tdm)
    error('To include photoselection weights, Sys.tdm must be given.');
  end
  if ischar(Exp.lightBeam)
    kLight = [0;1;0]; % beam propagating along yL
    switch Exp.lightBeam
      case 'perpendicular'
        alphaLight = -pi/2; % gives E-field along xL
      case 'parallel'
        alphaLight = pi; % gives E-field along zL
      case 'unpolarized'
        alphaLight = NaN; % unpolarized beam
      otherwise
        error('Unknown string in Exp.lightBeam. Use '''', ''perpendicular'', ''parallel'' or ''unpolarized''.');
    end
  else
    if ~iscell(Exp.lightBeam) || numel(Exp.lightBeam)~=2
      error('Exp.lightBeam should be a 2-element cell {k alpha}.')
    end
    kLight = Exp.lightBeam{1};  % propagation direction
    alphaLight = Exp.lightBeam{2};  % polarization angle
  end
end

% Process crystal orientations, crystal symmetry, and frame transforms
[Orientations,nOrientations,nSites,averageOverChi] = p_crystalorientations(Exp,Opt);


% Options
%---------------------------------------------------------------------
if ~isfield(Opt,'PerturbOrder'), Opt.PerturbOrder = 2; end

if (numel(Opt.PerturbOrder)~=1) || ~isreal(Opt.PerturbOrder)
  error('Opt.PerturbOrder must be either 1 or 2.');
end
switch Opt.PerturbOrder
  case 1, secondOrder = false;
  case 2, secondOrder = true;
  otherwise
    error('Opt.PerturbOrder must be either 1 or 2.');
end

if secondOrder
  logmsg(1,'2nd order perturbation theory');
else
  logmsg(1,'1st order perturbation theory');
end

if isfield(Opt,'ImmediateBinning') && Opt.ImmediateBinning
  error('Opt.ImmediateBinning is not supported for frequency-swept spectra.');
end
%---------------------------------------------------------------------


B0 = Exp.Field*1e-3; % mT -> T

gg = g*g.';
trgg = trace(gg);

% Prefactor for transition rate, one element per mS <-> mS-1 transition
mS_ = (S:-1:-S+1).';
c = bmagn/2 * sqrt(S*(S+1)-mS_.*(mS_-1));
c = c/planck/1e9;
c2 = c.^2;

nRows = nTransitions*nNucSublevels;
nu = zeros(nRows,nOrientations);
Intensity = zeros(nTransitions,nOrientations);
vecs = zeros(3,nOrientations);
E0 = zeros(1,nOrientations);

% Loop over all orientations
for iOri = 1:nOrientations
  R_L2M = erot(Orientations(iOri,:)).';  % lab frame -> molecular frame
  n0 = R_L2M*nB0_L;  % transform to molecular frame representation
  vecs(:,iOri) = n0;

  geff = norm(g.'*n0);
  E0_ = bmagn*geff*B0/planck/1e6; % MHz
  E0(iOri) = E0_;
  u = g.'*n0/geff; % molecular frame representation

  % Thermal polarization, using Zeeman level spacing
  if useTemperature
    % Levels ordered from mS = S down to mS = -S
    Populations = exp(-(2*S:-1:0).'*planck*E0_*1e6/boltzm/Exp.Temperature);
    Populations(isnan(Populations)) = 1; % T = 0: Inf*0 for ground state
    Populations = Populations/sum(Populations);
    Polarization = diff(Populations);
  else
    % no temperature: high-temperature limit, with kT replaced by h*nuRef/2
    nuRef = 1e3;  % reference frequency, MHz
    Polarization = E0_/(nuRef/2)/(2*S+1)*ones(nTransitions,1);
  end

  % Compute intensities
  %----------------------------------------------------------------

  % Compute photoselection weight if needed
  if usePhotoSelection
    if averageOverChi
      ori = Orientations(iOri,1:2);  % omit chi
    else
      ori = Orientations(iOri,1:3);
    end
    photoWeight = photoselect(Sys.tdm,ori,kLight,alphaLight);
    % Add isotropic contribution (from scattering)
    photoWeight = (1-Exp.lightScatter)*photoWeight + Exp.lightScatter;
  else
    photoWeight = 1;
  end

  % Compute quantum-mechanical transition rate
  if averageOverChi
    if mwmode.linearpolarizedMode
      TransitionRate = c2/2*(1-xi1^2)*(trgg-norm(g*u)^2);
    elseif mwmode.unpolarizedMode
      TransitionRate = c2/4*(1+xik^2)*(trgg-norm(g*u)^2);
    elseif mwmode.circpolarizedMode
      TransitionRate = c2/2*(1+xik^2)*(trgg-norm(g*u)^2) + ...
        mwmode.circSense*2*c2*xik*det(g)/geff;
    end
  else
    if mwmode.linearpolarizedMode
      nB1_ = R_L2M*nB1; % transform to molecular frame representation
      TransitionRate = c2*norm(cross(g.'*nB1_,u))^2;
    elseif mwmode.unpolarizedMode
      nk_ = R_L2M*nk; % transform to molecular frame representation
      TransitionRate = c2/2*(trgg-norm(g*u)^2-norm(cross(g.'*nk_,u))^2);
    elseif mwmode.circpolarizedMode
      nk_ = R_L2M*nk; % transform to molecular frame representation
      TransitionRate = c2*(trgg-norm(g*u)^2-norm(cross(g.'*nk_,u))^2) + ...
        mwmode.circSense*2*c2*xik*det(g)/geff;
    end
  end

  % Combine all factors into overall line intensity
  Intensity(:,iOri) = Polarization.*TransitionRate*photoWeight;

  % Compute line positions
  %----------------------------------------------------------------

  % Orientation-dependent quantities for zero-field splitting
  if highSpin
    Du = D*u;
    uDu = u.'*Du;
    uDDu = Du.'*Du;
    D1sq = uDDu - uDu^2;
    D2sq = 2*trDD + uDu^2 - 4*uDDu;
  end

  % Orientation-dependent quantities for hyperfine couplings
  if nNuclei>0
    nK = zeros(nNuclei,1);
    A1sq = zeros(nNuclei,1);
    A2 = zeros(nNuclei,1);
    A3 = zeros(nNuclei,1);
    DA = zeros(nNuclei,1);
    for n = 1:nNuclei
      K = A{n}*u;
      nK(n) = norm(K);
      k = K/nK(n);
      Ak = A{n}.'*k;
      kAu = Ak.'*u;
      kAAk = Ak.'*Ak;
      A1sq(n) = kAAk - kAu^2;
      A2(n) = detA(n)*(u.'*invA{n}*k);
      A3(n) = trAA(n) - nK(n)^2 - kAAk + kAu^2;
      if highSpin
        DA(n) = Du.'*Ak - uDu*kAu;
      end
    end
  end

  % Loop over all mS <-> mS-1 transitions
  for imS = 1:nTransitions
    mS = S + 1 - imS;

    % zeroth order
    dE = E0_;

    % first order
    if highSpin
      dE = dE - uDu/2*(3-6*mS);
    end
    if nNuclei>0
      dE = dE + mIc*nK;
    end

    % second order
    if secondOrder
      if highSpin
        x = D1sq*(4*S*(S+1)-3*(8*mS^2-8*mS+3))...
          - D2sq/4*(2*S*(S+1)-3*(2*mS^2-2*mS+1));
        dE = dE - x/(2*E0_);
      end
      if nNuclei>0
        x = mIc.^2*A1sq - (1-2*mS)*mIc*A2 + (II1-mIc.^2)*A3/2;
        dE = dE + x/(2*E0_);
        if highSpin
          y = (3-6*mS)*mIc*DA;
          dE = dE - y/E0_;
        end
      end
    end

    nu((imS-1)*nNucSublevels+(1:nNucSublevels),iOri) = dE;

  end

end

% Intensities
%-------------------------------------------------------------------
Int = repelem(Intensity,nNucSublevels,1)/nNucSublevels;

% Widths
%-------------------------------------------------------------------
Wid2 = zeros(nRows,nOrientations);  % squared FWHM, MHz^2

% H strain
if any(Sys.HStrain)
  Wid2 = Wid2 + Sys.HStrain.^2*vecs.^2;
end

% g strain and A strain (first nucleus only)
usegStrain = any(Sys.gStrain(:));
useAStrain = any(Sys.AStrain(:)) && nNuclei>0;
if usegStrain || useAStrain

  % g strain matrix (relative), scaled by E0 below
  if usegStrain
    gStrainMatrix = diag(Sys.gStrain(1,:)./Sys.g(1,:));
    if any(Sys.gFrame(1,:))
      R_g2M = erot(Sys.gFrame(1,:)).';  % g frame -> molecular frame
      gStrainMatrix = R_g2M*gStrainMatrix*R_g2M.';
    end
  else
    gStrainMatrix = zeros(3);
  end

  % Width^2 = n.'*(E0*G + c*mI*A)^2*n, with G = gStrainMatrix, A = AStrainMatrix
  qGG = sum(vecs.*(gStrainMatrix^2*vecs),1);
  lw2_g = E0.^2.*qGG;
  if useAStrain
    AStrainMatrix = diag(Sys.AStrain(1,:));  % MHz
    if any(Sys.AFrame(1,:))
      R_A2M = erot(Sys.AFrame(1,:)).'; % A frame -> molecular frame
      AStrainMatrix = R_A2M*AStrainMatrix*R_A2M.';
    end
    GA = gStrainMatrix*AStrainMatrix;
    qGA = sum(vecs.*((GA+GA.')*vecs),1);
    qAA = sum(vecs.*(AStrainMatrix^2*vecs),1);
    cmI = Sys.gAStrainCorr*mI{1}(:);
    lw2_gA = lw2_g + cmI*(E0.*qGA) + cmI.^2*qAA;  % one row per mI of first nucleus
    % mS outer, nuclear sublevels inner; first nucleus varies slowest
    n1 = numel(mI{1});
    idx = repmat(repelem(1:n1,nNucSublevels/n1),1,nTransitions);
    Wid2 = Wid2 + lw2_gA(idx,:);
  else
    Wid2 = Wid2 + lw2_g;
  end

end

% D strain
if any(Sys.DStrain(:))
  x = vecs(1,:);
  y = vecs(2,:);
  z = vecs(3,:);
  mS_ = (S:-1:-S).';
  mSS = mS_.^2-S*(S+1)/3;
  % Calculate derivatives of energy w.r.t. D and E
  dHdD_ = mSS*((3*z.^2-1)/2);
  dHdE_ = mSS*(3/2*(x.^2-y.^2));

  % Compute energy derivatives, pre-multiply with strain FWHMs.
  DeltaD = Sys.DStrain(1);
  DeltaE = Sys.DStrain(2);
  rDE = Sys.DStrainCorr; % correlation coefficient between D and E
  if rDE~=0
    % Transform correlated D-E strain to uncorrelated coordinates
    % Construct and diagonalize D-E covariance matrix
    R12 = rDE*DeltaD*DeltaE;
    CovMatrix = [DeltaD^2 R12; R12 DeltaE^2];
    [V,L] = eig(CovMatrix);
    L = sqrt(diag(L));
    dH1 = L(1)*(V(1,1)*dHdD_ + V(2,1)*dHdE_);
    dH2 = L(2)*(V(1,2)*dHdD_ + V(2,2)*dHdE_);
  else
    dH1 = dHdD_*DeltaD;
    dH2 = dHdE_*DeltaE;
  end

  % Calculate freq-domain linewidths, one row per mS <-> mS-1 transition
  lw1 = diff(dH1,1,1);
  lw2 = diff(dH2,1,1);
  Wid2 = Wid2 + repelem(lw1.^2+lw2.^2,nNucSublevels,1);
end

if any(Wid2(:))
  Wid = sqrt(Wid2);
else
  Wid = [];
end

% Transitions
%-------------------------------------------------------------------
% Levels are numbered by increasing energy; the first transition block
% (mS = S <-> S-1) involves the highest levels. Nuclear sublevel ordering
% is only approximate.
Transitions = zeros(nRows,2);
for imS = 1:nTransitions
  lowerLevels = (nTransitions-imS)*nNucSublevels + (1:nNucSublevels).';
  Transitions((imS-1)*nNucSublevels+(1:nNucSublevels),:) = ...
    [lowerLevels lowerLevels+nNucSublevels];
end

spec = 0;

% Reshape arrays in the case of crystals with site splitting
if nSites>1 && Opt.separateSites
  siz = [nRows*nSites, numel(nu)/nRows/nSites];
  nu = reshape(nu,siz);
  Int = reshape(Int,siz);
  if ~isempty(Wid), Wid = reshape(Wid,siz); end
end

% Arrange output
%---------------------------------------------------------------
Output = {nu,Int,Wid,Transitions,spec};
varargout = Output(1:max(nargout,1));

end
