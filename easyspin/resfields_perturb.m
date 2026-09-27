% resfields_perturb  Compute resonance fields for cw EPR 
%
%   ... = resfields_perturb(Sys,Exp)
%   ... = resfields_perturb(Sys,Exp,Opt)
%   [Pos,Int] = resfields_perturb(...)
%   [Pos,Int,Wid] = resfields_perturb(...)
%   [Pos,Int,Wid,Trans] = resfields_perturb(...)
%
%   Computes cw EPR line positions, intensities and widths using
%   perturbation theory.
%
%   Input:
%    Sys: spin system structure
%    Exp: experimental parameters
%      mwFreq              microwave frequency, in GHz
%      Range               sweep range, [sweepmin sweepmax], in mT
%      CenterSweep         sweep range, [center sweep], in mT
%                            negative fields are possible (field along -z(Lab))
%      Temperature         temperature, in K; if omitted: high-temperature limit
%      SampleFrame         Nx3 array of Euler angles (in radians) for sample/crystal orientations
%      CrystalSymmetry     crystal symmetry (space group etc.)
%      MolFrame            Euler angles (in radians) for molecular frame orientation
%      mwMode              excitation mode: 'perpendicular', 'parallel', {k pol}
%                            pol: polarization angle, 'circular+', 'circular-', 'unpolarized'
%    Opt: additional computational options
%      Verbosity           level of detail of printing; 0, 1, 2
%      PerturbOrder        perturbation order; 1 or 2
%      Sites               list of crystal sites to include (default []: all)
%
%   Output:
%    Pos     line positions (in mT)
%    Int     line intensities
%    Wid     Gaussian line widths, full width half maximum (FWHM), in mT
%    Trans   list of transitions (level indices; nuclear sublevel indices are approximate)

function varargout = resfields_perturb(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

% Compute resonance fields based on formulas from
% M. Iwasaki, J.Magn.Reson. 16, 417-423 (1974)
% Second-order perturbation treatment of the general spin hamiltonian in an
% arbitrary coordinate system
% https://doi.org/10.1016/0022-2364(74)90223-6

% Assert correct Matlab version
warning(chkmlver);

% Check number of input arguments.
switch nargin
  case 2, Opt = struct;
  case 3
  otherwise
    error('Use two or three inputs: resfields_perturb(Sys,Exp) or resfields_perturb(Sys,Exp,Opt)!');
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

if Sys.nElectrons~=1 
  err = sprintf('Perturbation theory available only for systems with 1 electron. Yours has %d.',Sys.nElectrons);
end
if any(Sys.L(:))
    err = sprintf('Perturbation theory not available for electron spin combined with orbital angular momentum!');
end
if any(Sys.AStrain)
%  err = ('A strain (Sys.AStrain) not supported with perturbation theory. Use matrix diagonalization or remove Sys.AStrain.');
end
if highSpin && any(Sys.DStrain(:)) && any(mod(Sys.S,1))
  err = ('D strain not supported for half-integer spins with perturbation theory. Use matrix diagonalization or remove Sys.DStrain.');
end
if any(Sys.DStrain(:)) && any(Sys.DFrame(:))
  err = 'D strain cannot be used with tilted D tensors.';
end
if isfield(Sys,'nn') && any(Sys.nn(:)~=0)
  err = 'Perturbation theory not available for nuclear-nuclear couplings (Sys.nn).';
end
if ~isempty(Sys.initState)
  err = 'Sys.initState is not supported by resfields_perturb.';
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

for iNuc = nNuclei:-1:1
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
DefaultExp.mwFreq = NaN;
DefaultExp.Range = NaN;
DefaultExp.CenterSweep = NaN;
DefaultExp.Temperature = NaN;
DefaultExp.mwMode = '';

DefaultExp.SampleFrame = [0 0 0];
DefaultExp.CrystalSymmetry = 1;
DefaultExp.MolFrame = [0 0 0];
DefaultExp.SampleRotation = [];

Exp = adddefaults(Exp,DefaultExp);

% Check for obsolete fields
if isfield(Exp,'CrystalOrientation')
  error('Exp.CrystalOrientation is no longer supported, use Exp.SampleFrame/Exp.MolFrame instead.');
end
if isfield(Exp,'Mode')
  error('Exp.Mode is no longer supported. Use Exp.mwMode instead.');
end

if isnan(Exp.mwFreq), error('Experiment.mwFreq is missing!'); end

% Sweep range from CenterSweep or Range (CenterSweep has precedence)
Exp.Range = p_sweeprange(Exp,false,true);
if isfield(Exp,'SearchRange'), Exp.Range = Exp.SearchRange; end
if isempty(Exp.Range), error('Exp.Range/Exp.CenterSweep is missing!'); end

% Negative fields: resonances at -B are obtained from those at +B via time
% reversal symmetry
mirror = any(Exp.Range<0);

% Determine excitation mode
[xi1,xik,nB1,nk,nB0_L,mwmode] = p_excitationgeometry(Exp.mwMode);

% Temperature, non-equilibrium populations
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

% Photoselection
if ~isfield(Exp,'lightBeam'), Exp.lightBeam = ''; end
if ~isfield(Exp,'lightScatter'), Exp.lightScatter = 0; end

usePhotoSelection = ~isempty(Exp.lightBeam) && Exp.lightScatter<1;

if usePhotoSelection
  if ~isfield(Sys,'tdm') || isempty(Sys.tdm)
    error('To include photoselection weights, Sys.tdm must be given.');
  end
  if ischar(Exp.lightBeam)
    k = [0;1;0]; % beam propagating along yL
    switch Exp.lightBeam
      case 'perpendicular'
        alpha = -pi/2; % gives E-field along xL
      case 'parallel'
        alpha = pi; % gives E-field along zL
      case 'unpolarized'
        alpha = NaN; % unpolarized beam
      otherwise
        error('Unknown string in Exp.lightBeam. Use '''', ''perpendicular'', ''parallel'' or ''unpolarized''.');
    end
    Exp.lightBeam = {k alpha};
  else
    if ~iscell(Exp.lightBeam) || numel(Exp.lightBeam)~=2
      error('Exp.lightBeam should be a 2-element cell {k alpha}.')
    end
  end
end

% Process crystal orientations, crystal symmetry, and frame transforms
[Orientations,nOrientations,nSites,averageOverChi] = p_crystalorientations(Exp,Opt);


% Options
%---------------------------------------------------------------------
if ~isfield(Opt,'Sites'), Opt.Sites = []; end
if ~isfield(Opt,'separateSites'), Opt.separateSites = true; end % internal; false keeps sites as columns

if ~isfield(Opt,'PerturbOrder'), Opt.PerturbOrder = 2; end

if numel(Opt.PerturbOrder)~=1 || ~isreal(Opt.PerturbOrder)
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


E0 = Exp.mwFreq*1e3; % MHz

for iNuc = nNuclei:-1:1
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
  E1A = zeros(max(nNucStates),nNuclei);
  Baxis = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);
  dB = Baxis(2)-Baxis(1);
  spec = zeros(1,Exp.nPoints);
end

if ~isnan(Exp.Temperature)
  Populations = exp(-planck*(2*S:-1:0).'*Exp.mwFreq*1e9/boltzm/Exp.Temperature);
  Populations = Populations/sum(Populations);
  Polarization = diff(Populations);
  Polarization = Polarization(end:-1:1);
else
  % no temperature: high-temperature limit, with kT replaced by h*nuRef/2
  nuRef = Exp.mwFreq;  % reference frequency, GHz
  Etrans = Exp.mwFreq;  % transition energy; all transitions are at mwFreq, GHz
  Polarization = Etrans/(nuRef/2)/(2*S+1)*ones(2*S,1);
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
  R_L2M = erot(Orientations(iOri,:)).';  % lab frame -> molecular frame
  n0 = R_L2M*nB0_L;  % transform to molecular frame representation
  vecs(:,iOri) = n0;
  
  geff(iOri) = norm(g.'*n0);
  u = g.'*n0/geff(iOri); % molecular frame representation
  
  % frequency to field conversion factor
  preOri = 1e6*planck/(geff(iOri)*bmagn);
  
  % Compute intensities
  %----------------------------------------------------------------
  % Compute photoselection weight if needed
  if usePhotoSelection
    k = Exp.lightBeam{1};  % propagation direction
    alpha = Exp.lightBeam{2};  % polarization angle
    if averageOverChi
      ori = Orientations(iOri,1:2);  % omit chi
    else
      ori = Orientations(iOri,1:3);
    end
    photoWeight = photoselect(Sys.tdm,ori,k,alpha);
    % Add isotropic contribution (from scattering)
    photoWeight = (1-Exp.lightScatter)*photoWeight + Exp.lightScatter;
  else
    photoWeight = 1;
  end

  % Compute quantum-mechanical transition rate
  if averageOverChi
    if mwmode.linearpolarizedMode
      TransitionRate(:,iOri) = c2/2*(1-xi1^2)*(trgg-norm(g*u)^2);
    elseif mwmode.unpolarizedMode
      TransitionRate(:,iOri) = c2/4*(1+xik^2)*(trgg-norm(g*u)^2);
    elseif mwmode.circpolarizedMode
      TransitionRate(:,iOri) = c2/2*(1+xik^2)*(trgg-norm(g*u)^2);
      circularTerm = mwmode.circSense*2*c2*xik*det(g)/norm(g.'*n0);
    end
  else
    if mwmode.linearpolarizedMode
      nB1_ = R_L2M*nB1; % transform to molecular frame representation
      TransitionRate(:,iOri) = c2*norm(cross(g.'*nB1_,u))^2;
    elseif mwmode.unpolarizedMode
      nk_ = R_L2M*nk; % transform to molecular frame representation
      TransitionRate(:,iOri) = c2/2*(trgg-norm(g*u)^2-norm(cross(g.'*nk_,u))^2);
    elseif mwmode.circpolarizedMode
      nk_ = R_L2M*nk; % transform to molecular frame representation
      TransitionRate(:,iOri) = c2*(trgg-norm(g*u)^2-norm(cross(g.'*nk_,u))^2);
      circularTerm = mwmode.circSense*2*c2*det(g)*xik/norm(g.'*n0);
    end
  end
  if mwmode.circpolarizedMode
    % circular term changes sign for mirrored resonances at -B
    TransitionRateMirror = TransitionRate(:,iOri) - circularTerm;
    TransitionRate(:,iOri) = TransitionRate(:,iOri) + circularTerm;
  else
    TransitionRateMirror = TransitionRate(:,iOri);
  end

  % Compute Aasa-Vänngård 1/g factor (frequency-to-field conversion factor)
  dBdE = (planck/bmagn*1e9)/geff(iOri);
  
  % Combine all factors into overall line intensity
  Intensity(:,iOri) = Polarization.*TransitionRate(:,iOri)*dBdE*photoWeight;
  IntensityMirror(:,iOri) = Polarization.*TransitionRateMirror*dBdE*photoWeight;
  
  if highSpin
    Du = D*u;
    uDu = u.'*Du;
    uDDu = Du.'*Du;
    D1sq = uDDu - uDu^2;
    D2sq = 2*trDD + uDu^2 - 4*uDDu;
  end

  imS = 0;
  for mS = S:-1:-S+1
    imS = imS + 1;

    % first-order correction
    if highSpin
      E1D = -uDu/2*(3-6*mS);
    else
      E1D = 0;
    end
    if nNuclei>0
      for iNuc = 1:nNuclei
        K = A{iNuc}*u;
        nK = norm(K);
        k(:,iNuc) = K/nK;
        E1A_ = mI{iNuc}*nK;
        if immediateBinning
          E1A(1:nNucStates(iNuc),iNuc) = E1A_(:);
        else
          E1A(:,iNuc) = E1A_(idxn(:,iNuc)).';
        end
      end
    else
      E1A = 0;
    end

    % second-order correction
    E2D = 0;
    if secondOrder
      if highSpin
        x =  D1sq*(4*S*(S+1)-3*(8*mS^2-8*mS+3))...
          - D2sq/4*(2*S*(S+1)-3*(2*mS^2-2*mS+1));
        E2D = -x./(2*E0);
      else
        E2DA = 0;
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
          E2A_ = +x./(2*E0);
          if immediateBinning
            E2A(1:nNucStates(n),n) = E2A_(:);
          else
            E2A(:,n) = E2A_(idxn(:,n)).';
          end
          if highSpin
            DA = Du.'*Ak - uDu*kAu;
            y = DA*(3-6*mS)*mI{n};
            E2DA_ = -y./E0;
            if immediateBinning
              E2DA(1:nNucStates(n),n) = E2DA_;
            else
              E2DA(:,n) = E2DA_(idxn(:,n)).';
            end
          else
            E2DA = 0;
          end
        end
      else
        E2DA = 0;
        E2A = 0;
      end
    else
      E2A = 0;
      E2DA = 0;
    end
    
    if immediateBinning
      B0 = (E0-E1D-E2D)*preOri*1e3; % mT
      % compute B shifts
      Bshifts = (-(E1A+E2A+E2DA))*preOri*1e3; % mT
      % directly accumulate into spectrum
      spec = spec + Intensity(imS,iOri)*Exp.AccumWeights(iOri)*...
        multinucstick(B0,nNucStates,Bshifts,Baxis(1),dB,Exp.nPoints);
      if mirror
        spec = spec + IntensityMirror(imS,iOri)*Exp.AccumWeights(iOri)*...
          multinucstick(-B0,nNucStates,-Bshifts,Baxis(1),dB,Exp.nPoints);
      end
    else
      if secondOrder
        Bfinal{imS}(iOri,:) = (E0-E1D-E2D-sum(E1A+E2A+E2DA,2))*preOri;
      else
        Bfinal{imS}(iOri,:) = (E0-E1D-sum(E1A,2))*preOri;
      end
    end
    
  end
  
end

if immediateBinning
  B = [];
  Int = [];
  Wid = [];
  Transitions = [];
  spec = spec/dB/prod(nNucStates);
  spec = spec*(2*pi); % powder chi integral
else
  % Positions
  %-------------------------------------------------------------------
  B = [Bfinal{:}].';
  B = B*1e3;  % T -> mT
  
  % Intensities
  %-------------------------------------------------------------------
  nNucSublevels = prod(nNucStates);
  Int = repelem(Intensity,nNucSublevels,1)/nNucSublevels;
  Int = flipud(Int);
  
  % Widths
  %-------------------------------------------------------------------
  if any(Sys.HStrain)
    lw2 = sum(Sys.HStrain.^2*vecs.^2,1); % MHz^2
    lw = sqrt(lw2)*1e6*planck./geff/bmagn*1e3; % mT
    Wid = repmat(lw,nNucTrans*2*S,1);
  else
    Wid = 0;
  end
  
  if any(Sys.gStrain(:)) || any(Sys.AStrain(:))
    
    if any(Sys.gStrain(:))
      gStrainMatrix = diag(Sys.gStrain(1,:)./Sys.g(1,:))*Exp.mwFreq*1e3;  % -> MHz
      if any(Sys.gFrame(:))
        R_g2M = erot(Sys.gFrame(1,:)).';  % g frame -> molecular frame
        gStrainMatrix = R_g2M*gStrainMatrix*R_g2M.';
      end
    else
      gStrainMatrix = zeros(3);
    end
    
    if any(Sys.AStrain) && Sys.nNuclei>0
      AStrainMatrix = diag(Sys.AStrain);
      if any(Sys.AFrame(1,:))
        Rp = erot(Sys.AFrame(1,:)).'; % A frame -> molecular frame
        AStrainMatrix = Rp*AStrainMatrix*Rp.';
      end
      corr = Sys.gAStrainCorr;
      mI1 = -Sys.I(1):+Sys.I(1);
      lw2 = zeros(numel(mI1),nOrientations);
      for idx = 1:numel(mI1)
        StrainMatrix = gStrainMatrix + corr*(mI1(idx))*AStrainMatrix;
        lw2(idx,:) = sum(vecs.*(StrainMatrix^2*vecs),1);
      end
      Wid_gA = sqrt(lw2)*planck*1e6./repmat(geff,numel(mI1),1)/bmagn*1e3; % MHz -> mT
      % mS outer, nuclear sublevels inner; first nucleus varies slowest
      idx = repmat(repelem(1:numel(mI1),nNucTrans/numel(mI1)),1,2*S);
      Wid = sqrt(Wid_gA(idx,:).^2 + Wid.^2);
    else
      lw2 = sum(vecs.*(gStrainMatrix^2*vecs),1);
      Wid_gA = sqrt(lw2)*planck*1e6./geff/bmagn*1e3; % MHz -> mT
      Wid = sqrt(repmat(Wid_gA.^2,2*S*nNucTrans,1)+Wid.^2);
    end

  end

  if any(Sys.DStrain(:))
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
      dH1 = L(1)*(V(1,1)*dHdD_ + V(2,1)*dHdE_);
      dH2 = L(2)*(V(1,2)*dHdD_ + V(2,2)*dHdE_);
    else
      dH1 = dHdD_*DeltaD;
      dH2 = dHdE_*DeltaE;
    end

    % Calculate freq-domain linewidths, one row per mS <-> mS-1 transition
    lw1 = diff(dH1,1,1);
    lw2 = diff(dH2,1,1);
    % convert from MHz to mT
    MHz2mT = (planck/bmagn*1e9)./geff;
    lw1 = lw1.*MHz2mT;
    lw2 = lw2.*MHz2mT;
    Wid2_DE = repelem(lw1.^2+lw2.^2,nNucTrans,1);
    Wid = sqrt(Wid2_DE + Wid.^2);
  end

  if ~any(Wid(:))
    Wid = [];
  end

  % Transitions
  %-------------------------------------------------------------------
  % Levels are numbered by increasing energy; the first transition block
  % (mS = S <-> S-1) involves the highest levels. Nuclear sublevel ordering
  % is only approximate.
  Transitions = zeros(2*S*nNucSublevels,2);
  for imS = 1:2*S
    lowerLevels = (2*S-imS)*nNucSublevels + (1:nNucSublevels).';
    Transitions((imS-1)*nNucSublevels+(1:nNucSublevels),:) = ...
      [lowerLevels lowerLevels+nNucSublevels];
  end

  % Negative fields: add mirrored resonances at -B, and remove resonances
  % without any position in range
  if mirror
    IntMirror = repelem(IntensityMirror,nNucSublevels,1)/nNucSublevels;
    IntMirror = flipud(IntMirror);
    B = [B; -B];
    Int = [Int; IntMirror];
    if ~isempty(Wid), Wid = [Wid; Wid]; end
    Transitions = [Transitions; Transitions];
    keep = any(B>=Exp.Range(1) & B<=Exp.Range(2),2);
    B = B(keep,:);
    Int = Int(keep,:);
    if ~isempty(Wid), Wid = Wid(keep,:); end
    Transitions = Transitions(keep,:);
  end

  spec = 0;

end

% Reshape arrays in the case of crystals with site splitting
if nSites>1 && Opt.separateSites
  nRows = size(B,1);
  siz = [nRows*nSites, numel(B)/nRows/nSites];
  B = reshape(B,siz);
  if ~isempty(Int), Int = reshape(Int,siz); end
  if ~isempty(Wid), Wid = reshape(Wid,siz); end
end

% Arrange output
%---------------------------------------------------------------
Output = {B,Int,Wid,Transitions,spec};
varargout = Output(1:max(nargout,1));

end
