% resfreqs_perturb Compute resonance frequencies for frequency-swept EPR
%
%   ... = resfreqs_perturb(Sys,Exp)
%   ... = resfreqs_perturb(Sys,Exp,Opt)
%   [Pos,Int] = resfreqs_perturb(...)
%   [Pos,Int,Wid] = resfreqs_perturb(...)
%   [Pos,Int,Wid,Trans] = resfreqs_perturb(...)
%
%   Computes frequency-domain EPR line positions, intensities and widths using
%   perturbation theory.
%
%   Input:
%    Sys: spin system structure
%    Exp: experimental parameter settings
%      Field               static field, in mT
%      Range               frequency sweep range, [numin numax], in GHz
%      CenterField         frequency sweep range, [center sweep], in GHz
%      Temperature         temperature, in K
%      CrystalOrientation  nx3 array of Euler angles (in radians) for crystal orientations
%      CrystalSymmetry     crystal symmetry (space group etc.)
%      MolFrame            Euler angles (in radians) for molecular frame orientation
%      mwPolarization      'linear', 'circular+', 'circular-', 'unpolarized'
%      Mode                excitation mode: 'perpendicular', 'parallel', [k_tilt alpha_pol]
%    Opt: additional computational options
%      Verbosity           level of detail of printing; 0, 1, 2
%      PerturbOrder        perturbation order; 1 or 2
%      Sites               list of crystal sites to include (default []: all)
%
%   Output:
%    Pos     line positions (in mT)
%    Int     line intensities
%    Wid     Gaussian line widths, full width half maximum (FWHM)
%    Trans   list of transitions included in the computation

function varargout = resfreqs_perturb(Sys,Exp,Opt)

% Compute resonance fields based on formulas from Iwasaki, J.Magn.Reson. 16, 417-423 (1974)

% Assert correct Matlab version
error(chkmlver);

% Check number of input arguments.
switch nargin
  case 0, help(mfilename); return;
  case 2, Opt = struct;
  case 3,
  otherwise
    error('Use two or three inputs: refields_perturb(Sys,Exp) or refields_perturb(Sys,Exp,Opt)!');
end

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;

% Spin system
%------------------------------------------------------
[Sys,err] = validatespinsys(Sys);
error(err);
S = Sys.S;
highSpin = any(S>1/2);

if Sys.nElectrons~=1
  err = sprintf('Perturbation theory available only for systems with 1 electron. Yours has %d.',Sys.nElectrons);
end
if any(Sys.L(:))
  err = sprintf('Pertrubation theory not available for electron spin coupled to orbital angular momentum!');
end
if any(Sys.AStrain)
%  err = ('A strain (Sys.AStrain) not supported with perturbation theory. Use matrix diagonalization or remove Sys.AStrain.');
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
if isfield(Sys,'Pop') && any(Sys.Pop(:))
  err = 'Sys.Pop is not supported by resfreqs_perturb.';
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
  if ~Sys.fullD
    R_D2M = erot(Sys.DFrame).'; % D frame -> molecular frame
    D = R_D2M*diag(Sys.D)*R_D2M.';
  end
  % make D traceless (required for Iwasaki expressions)
  D = D - eye(3)*trace(D)/3;
end

nTransitions = 2*S; % number of allowed transitions for one electron spin

I = Sys.I;
nNuclei = Sys.nNuclei;
if nNuclei>0
  nNucStates = 2*I+1;
else
  nNucStates = 1;
end

% Guard against zero hyperfine couplings
% (otherwise inv(A) gives error further down)
if nNuclei>0
  if ~Sys.fullA
    if any(Sys.A(:)==0)
      error('All hyperfine coupling constants must be non-zero.');
    end
  end
end

for iNuc = 1:nNuclei
  if Sys.fullA
    A{iNuc} = Sys.A((iNuc-1)*3+(1:3),:);
  else
    R_A2M = erot(Sys.AFrame(iNuc,:)).'; % A frame -> molecular frame
    A_ = diag(Sys.A(iNuc,:));
    A{iNuc} = R_A2M*A_*R_A2M.';
  end
  mI{iNuc} = -I(iNuc):I(iNuc);
  idxn{iNuc} = 1:nNucStates(iNuc);
end


% Experiment
%------------------------------------------------------
DefaultExp.Range = NaN;
DefaultExp.Field = NaN;
DefaultExp.CenterSweep = NaN;
DefaultExp.Temperature = NaN;
DefaultExp.Mode = 'perpendicular';
DefaultExp.mwPolarization = '';

DefaultExp.CrystalOrientation = [];
DefaultExp.CrystalSymmetry = '';
DefaultExp.MolFrame = [];

Exp = adddefaults(Exp,DefaultExp);

err = '';
if ~isfield(Exp,'Field'), err = 'Exp.Field is missing.'; end

p_excitationgeometry;

if isfield(Exp,'Temperature')
  if numel(Exp.Temperature)>1
    err = 'Exp.Temperature must be a single number.';
  end
  if isinf(Exp.Temperature)
    err = 'If given, Exp.Temperature must have a finite value.';
  end
else
  Exp.Temperature = NaN;
end
error(err);

if ~isfield(Opt,'Sites')
  Opt.Sites = [];
end


% Process crystal orientations, crystal symmetry, and frame transforms
[Orientations,nOrientations,nSites,AverageOverChi] = p_crystalorientations(Exp,Opt);


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

if ~isfield(Opt,'ImmediateBinning'), Opt.ImmediateBinning = 0; end
%---------------------------------------------------------------------


B0 = Exp.Field*1e-3; % mT -> T

for iNuc = 1:nNuclei
  A_ = A{iNuc};
  detA(iNuc) = det(A_);
  invA{iNuc} = inv(A_); % gives an error with zero hf couplings
  trAA(iNuc) = trace(A_.'*A_);
end

if highSpin
  trDD = trace(D^2);
end
II1 = I.*(I+1);

immediateBinning = Opt.ImmediateBinning;

if immediateBinning
else
  if nNuclei>0
    idxn = allcombinations(idxn{:});
  else
    idxn = 1;
  end
  nNucTrans = size(idxn,1);
end

if immediateBinning
  dE1A = zeros(max(nNucStates),nNuclei);
  nuaxis = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);
  dnu = nuaxis(2)-nuaxis(1);
  spec = zeros(1,Exp.nPoints);
end

if ~isnan(Exp.Temperature)
  Populations = exp(-planck*(2*S:-1:0).'*Exp.mwFreq*1e9/boltzm/Exp.Temperature);
  Populations = Populations/sum(Populations);
  Polarization = diff(Populations);
  Polarization = Polarization(end:-1:1);
else
  Polarization = ones(2*S,1);
end

gg = g*g.';
trgg = trace(gg);

% prefactor for transition rate
mS = S:-1:-S+1;
c = bmagn/2 * sqrt(S*(S+1)-mS.*(mS-1));
c = c/planck/1e9;
c2 = c.^2;

% Loop over all orientations
for iOri = nOrientations:-1:1
  R = erot(Orientations(iOri,:));
  n0 = R.'*nB0;  % transform to molecular frame representation
  vecs(:,iOri) = n0;
  
  geff(iOri) = norm(g.'*n0);
  E0 = bmagn*geff(iOri)*B0/planck/1e6; % MHz
  u = g.'*n0/geff(iOri); % molecular frame representation
    
  % Compute intensities
  %----------------------------------------------------------------

  % Compute quantum-mechanical transition rate
  if AverageOverChi
    if linearpolarizedMode
      TransitionRate(:,iOri) = c2/2*(1-xi1^2)*(trgg-norm(g*u)^2);
    elseif unpolarizedMode
      TransitionRate(:,iOri) = c2/4*(1+xik^2)*(trgg-norm(g*u)^2);
    elseif circpolarizedMode
      TransitionRate(:,iOri) = c2/2*(1+xik^2)*(trgg-norm(g*u)^2) + ...
        circSense*2*c2*xik^2*det(g)/norm(g.'*n0);
    end
  else
    if linearpolarizedMode
      nB1_ = R.'*nB1; % transform to molecular frame representation
      TransitionRate(:,iOri) = c2*norm(cross(g.'*nB1_,u))^2;
    elseif unpolarizedMode
      nk_ = R.'*nk; % transform to molecular frame representation
      TransitionRate(:,iOri) = c2/2*(trgg-norm(g*u)^2-norm(cross(g.'*nk_,u))^2);
    elseif circpolarizedMode
      nk_ = R.'*nk; % transform to molecular frame representation
      TransitionRate(:,iOri) = c2*(trgg-norm(g*u)^2-norm(cross(g.'*nk_,u))^2) + ...
        circSense*2*c2*det(g)*xik/norm(g.'*n0);
    end
  end

  % Combine all factors into overall line intensity
  Intensity(:,iOri) = Polarization.*TransitionRate(:,iOri);
  
  if highSpin
    Du = D*u;
    uDu = u.'*Du;
    uDDu = Du.'*Du;
    D1sq = uDDu - uDu^2;
    D2sq = 2*trDD + uDu^2 - 4*uDDu;
  end

  imS = 0;
  for mS = S:-1:-S+1
    % for mS <-> mS-1 transition
    imS = imS + 1;

    % first-order
    if highSpin
      dE1D = -uDu/2*(3-6*mS);
    else
      dE1D = 0;
    end
    if nNuclei>0
      for iNuc = 1:nNuclei
        K = A{iNuc}*u;
        nK = norm(K);
        k(:,iNuc) = K/nK;
        dE1A_ = mI{iNuc}*nK;
        if immediateBinning
          dE1A(1:nNucStates(iNuc),iNuc) = dE1A_(:);
        else
          dE1A(:,iNuc) = dE1A_(idxn(:,iNuc)).';
        end
      end
    else
      dE1A = 0;
    end

    % second-order
    dE2D = 0;
    if secondOrder
      if highSpin
        x =  D1sq*(4*S*(S+1)-3*(8*mS^2-8*mS+3))...
          - D2sq/4*(2*S*(S+1)-3*(2*mS^2-2*mS+1));
        dE2D = -x./(2*E0);
      else
        dE2DA = 0;
      end
      if nNuclei>0
        for n = 1:nNuclei
          k_ = k(:,n);
          Ak = A{n}.'*k_;
          kAu = Ak.'*u;
          kAAk = norm(Ak)^2;
          A1sq = kAAk - kAu^2;
          A2 = detA(n)*(u.'*invA{n}*k_);
          A3 = trAA(n) - norm(A{n}*u)^2 - kAAk + kAu^2;
          x = A1sq*mI{n}.^2 - A2*(1-2*mS)*mI{n} + A3/2*(II1(n)-mI{n}.^2);
          dE2A_ = +x./(2*E0);
          if immediateBinning
            dE2A(1:nNucStates(n),n) = dE2A_(:);
          else
            dE2A(:,n) = dE2A_(idxn(:,n)).';
          end
          if highSpin
            DA = Du.'*Ak - uDu*kAu;
            y = DA*(3-6*mS)*mI{n};
            dE2DA_ = -y./E0;
            if immediateBinning
              dE2DA(1:nNucStates(n),n) = dE2DA_;
            else
              dE2DA(:,n) = dE2DA_(idxn(:,n)).';
            end
          else
            dE2DA = 0;
          end
        end
      else
        dE2DA = 0;
        dE2A = 0;
      end
    else
      dE2A = 0;
      dE2DA = 0;
    end
    
    if immediateBinning
      dE = dE0 + dE1D + dE2D + dE1A + dE2A + dE2DA;
      % directly accumulate into spectrum
      spec = spec + Intensity(imS,iOri)*Exp.AccumWeights(iOri)*...
        multinucstick(dE,nNucStates,Bshifts,nuaxis(1),dnu,Exp.nPoints);
    else
      if secondOrder
        dEfinal{imS}(iOri,:) = E0+dE1D+dE2D+sum(dE1A+dE2A+dE2DA,2);
      else
        dEfinal{imS}(iOri,:) = E0+dE1D+sum(dE1A,2);
      end
    end
    
  end
  
end

if immediateBinning
  nu = [];
  Int = [];
  Wid = [];
  Transitions = [];
  spec = spec/dnu/prod(nNucStates);
  spec = spec*(2*pi); % powder chi integral
else
  % Positions
  %-------------------------------------------------------------------
  nu = [];
  for iTrans = 1:nTransitions
    nu = [nu dEfinal{iTrans}];
  end
  nu = nu.';
  
  % Intensities
  %-------------------------------------------------------------------
  nNucSublevels = prod(nNucStates);
  Int = [];
  for iTrans = nTransitions:-1:1
    Int = [Int; repmat(Intensity(iTrans,:),nNucSublevels,1)];
  end
  Int = Int/nNucSublevels;
  
  % Widths
  %-------------------------------------------------------------------
  if any(Sys.HStrain)
    lw2 = sum(Sys.HStrain.^2*vecs.^2,1);
    lw = sqrt(lw2);
    Wid = repmat(lw,nNucTrans*2*S,1);
  else
    Wid = 0;
  end
    
  if any(Sys.gStrain(:)) || any(Sys.AStrain(:))
    
    if any(Sys.gStrain(:))
      gStrainMatrix = diag(Sys.gStrain(1,:)./Sys.g(1,:))*E0*1e3; % -> MHz
      if any(Sys.gFrame(:))
        R_g2M = erot(Sys.gFrame(1,:)).'; % g frame -> molecular frame
        gStrainMatrix = R_g2M*gStrainMatrix*R_g2M.';
      end
    else
      gStrainMatrix = zeros(3);
    end
    
    if any(Sys.AStrain(:)) && (Sys.nNuclei>0)
      AStrainMatrix = diag(Sys.AStrain);
      if isfield(Sys,'AFrame')
        R_A2M = erot(Sys.AFrame(1,:)).'; % A frame -> molecular frame
        AStrainMatrix = R_A2M*AStrainMatrix*R_A2M.';
      end
      corr = Sys.gAStrainCorr;
      mI1 = -Sys.I(1):+Sys.I(1);
      for idx = 1:numel(mI1)
        StrainMatrix = gStrainMatrix + corr*(mI1(idx))*AStrainMatrix;
        for iOri = 1:nOrientations
          lw2(idx,iOri) = vecs(:,iOri).'*StrainMatrix.^2*vecs(:,iOri);
        end
      end
      Wid_gA = sqrt(lw2); % MHz
      idx = repmat(1:numel(mI1),2*S*nNucTrans/numel(mI1),1);
      Wid = Wid_gA(idx(:),:);
    else
      StrainMatrix = gStrainMatrix;
      for iOri = 1:nOrientations
        lw2(1,iOri) = vecs(:,iOri).'*StrainMatrix.^2*vecs(:,iOri);
      end
      Wid_gA = sqrt(lw2); % MHz
      Wid = repmat(Wid_gA,2*S*nNucTrans,1);
    end
  
  elseif any(Sys.DStrain(:))
    if any(Sys.DFrame(:))
      error('Cannot use D/E strain with tilted D tensor.');
    end
    x = vecs(1,:);
    y = vecs(2,:);
    z = vecs(3,:);
    mS = (S:-1:-S).';
    mSS = mS.^2-S*(S+1)/3;
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
      dHdD_ = L(1)*(V(1,1)*dHdD_ + V(1,2)*dHdE_);
      dHdE_ = L(2)*(V(2,1)*dHdD_ + V(2,2)*dHdE_);
    else
      dHdD_ = dHdD_*DeltaD;
      dHdE_ = dHdE_*DeltaE;
    end
    
    % Calculate freq-domain linedwidths
    lwD = diff(dHdD_,1,1);
    lwE = diff(dHdE_,1,1);
    Wid2_DE = repmat(lwD.^2+lwE.^2,nNucTrans,1);
    Wid = sqrt(Wid2_DE + Wid.^2);
    
  else
    Wid = [];
  end
  
  % Transitions
  %-------------------------------------------------------------------
  nI = prod(nNucStates);
  Transitions = [];
  Manifold = (1:nI).';
  for k = 1:2*S
    Transitions = [Transitions; [Manifold Manifold+nI]];
    Manifold = Manifold + nI;
  end
  
  spec = 0;
end

% Reshape arrays in the case of crystals with site splitting
d = dbstack;
pepperCall = numel(d)>2 && strcmp(d(2).name,'pepper');
if (nSites>1) && ~pepperCall
  siz = [nTransitions*nSites, numel(nu)/nTransitions/nSites];
  nu = reshape(nu,siz);
  if ~isempty(Int), Int = reshape(Int,siz); end
  if ~isempty(Wid), Wid = reshape(Wid,siz); end
end

% Arrange output
%---------------------------------------------------------------
Output = {nu,Int,Wid,Transitions,spec};
varargout = Output(1:max(nargout,1));

return

