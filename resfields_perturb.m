function varargout = resfields_perturb(Sys,Exp,Opt)

% Compute resonance fields based on formulas from
% Iwasaki, J.Magn.Reson. 16, 417-423 (1974)

% Assert correct Matlab version
error(chkmlver);

% Check number of input arguments.
switch (nargin)
  case 0, help(mfilename); return;
  case 2, Opt = struct('unused',NaN);
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
if any(Sys.AStrain)
%  err = ('A strain (Sys.AStrain) not supported with perturbation theory. Use matrix diagonalization or remove Sys.AStrain.');
end
if highSpin && any(Sys.DStrain) && any(mod(Sys.S,1))
  err = ('D strain not supported for half-integer spins with perturbation theory. Use matrix diagonalization or remove Sys.DStrain.');
end
if any(Sys.DStrain(:)) && any(Sys.Dpa(:))
  err = 'D stain cannot be used with tilted D tensors.';
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
  Rg = erot(Sys.gpa);
  g = Rg*diag(Sys.g)*Rg.';
end

if highSpin
  if ~Sys.fullD
    RD = erot(Sys.Dpa);
    D = RD*diag(Sys.D)*RD.';
  end
  % make D traceless (required for Iwasaki expressions)
  D = D - eye(3)*trace(D)/3;
end

nTransitions = 2*S; % number of allowed transitions for one electron spin

I = Sys.I;
nNuclei = Sys.nNuclei;
if (nNuclei>0)
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
    RA = erot(Sys.Apa(iNuc,:));
    A_ = diag(Sys.A(iNuc,:))*Sys.Ascale(iNuc);
    A{iNuc} = RA*A_*RA';
  end
  mI{iNuc} = -I(iNuc):I(iNuc);
  idxn{iNuc} = 1:nNucStates(iNuc);
end


% Experiment
%------------------------------------------------------
if ~isfield(Exp,'CrystalSymmetry'), Exp.CrystalSymmetry = ''; end

err = '';
if ~isfield(Exp,'mwFreq'), err = 'Exp.mwFreq is missing.'; end
if ~isfield(Exp,'Orientations'), err = 'Exp.Orientations is missing'; end
if isfield(Exp,'Detection'), err = 'Exp.Detection is obsolete. Use Exp.Mode instead.'; end

if isfield(Exp,'Temperature')
  if numel(Exp.Temperature)~=1
    err = 'Exp.Temperature must be a single number.';
  end
  if isinf(Exp.Temperature)
    err = 'If given, Exp.Temperature must have a finite value.';
  end
else
  Exp.Temperature = NaN;
end
error(err);

% Orientations
%------------------------------------------------------
Orientations = Exp.Orientations;
[n1,n2] = size(Orientations);
if ((n2==2)||(n2==3)) && (n1~=2) && (n1~=3)
  Orientations = Orientations.';
end
[nAngles,nOrientations] = size(Orientations);
switch nAngles
  case 2
    AverageOverChi = true;
    Orientations(3,end) = 0; % Entire chi column is set to 0.
  case 3
    AverageOverChi = false;
  otherwise
    error('Orientations array has %d rows instead of 2 or 3.',nAngles);
end

% Add symmetry-related sites if space group symmetry is given
if ~isempty(Exp.CrystalSymmetry)
  R = sitetransforms(Exp.CrystalSymmetry);
  nSites  = numel(R);
  allOrientations = zeros(nOrientations*nSites,3);
  idx = 1;
  for iOri = 1:nOrientations
    xyz0 = erot(Orientations(:,iOri)).'; % xL, yL, zL along columns
    for iSite = 1:nSites
      xyz = R{iSite}*xyz0; % active rotation
      allOrientations(idx,:) = eulang(xyz.',1);
      idx = idx + 1;
    end
  end
  Orientations = allOrientations.';
  [nAngles,nOrientations] = size(Orientations);
else
  nSites = 1;
end


% Options
%----------------------------------------------------------
if ~isfield(Opt,'PerturbOrder'), Opt.PerturbOrder = 2; end
if ~isfield(Opt,'DirectAccumulation'), Opt.DirectAccumulation = 0; end
secondOrder = (Opt.PerturbOrder==2);
if secondOrder
  logmsg(1,'2nd order perturbation theory');
else
  logmsg(1,'1st order perturbation theory');
end
if (Opt.PerturbOrder>2)
  error('Only 1st and 2nd order perturbation theory are supported.');
end
%----------------------------------------------------------


E0 = Exp.mwFreq*1e3; % MHz

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

directAccumulation = Opt.DirectAccumulation;

if directAccumulation
else
  if (nNuclei>0)
    idxn = allcombinations(idxn{:});
  else
    idxn = 1;
  end
  nNucTrans = size(idxn,1);
end

if directAccumulation
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
  Polarization = ones(2*S,1);
end

gg = g*g.';
trgg = trace(gg);

% prefactor for transition rate
mS = S:-1:-S+1;
c = (bmagn/planck/1e9/2)^2 * (S*(S+1)-mS.*(mS-1));

% Loop over all orientations
for iOri = nOrientations:-1:1
  [n1x,n1y,n0] = erot(Orientations(:,iOri));
  n0 = n0(:);
  vecs(:,iOri) = n0;
  
  geff(iOri) = norm(g.'*n0);
  u = g.'*n0/geff(iOri);
  
  % frequency to field conversion factor
  preOri = 1e6*planck/(geff(iOri)*bmagn);
  
  % Compute intensities
  %----------------------------------------------------------------
  % Transition rate
  % - see Weil/Bolton p.104
  % - A.Lund et al, 2008, appendix, eq. (A5)
  % - Iwasaki 1974 eq. [40]
  % - Kneubuehl 1961 eq. [5]
  if ParallelMode
    TransitionRate(:,iOri) = 0;
  else
    if (AverageOverChi)
      TransitionRate(:,iOri) = c*(trgg-u.'*gg*u)/2; % Iwasaki
    else
      TransitionRate(:,iOri) = c*norm(cross(g.'*n1x.',u))^2; % Kneubuehl
    end
  end
  
  % Aasa-Vänngård 1/g factor
  dBdE = (planck/bmagn*1e9)/geff(iOri);
  
  % Combine all factors into overall line intensity
  Intensity(:,iOri) = Polarization.*TransitionRate(:,iOri)*dBdE;
  
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

    % first-order
    if highSpin
      E1D = -uDu/2*(3-6*mS);
    else
      E1D = 0;
    end
    if (nNuclei>0)
      for iNuc = 1:nNuclei
        K = A{iNuc}*u;
        nK = norm(K);
        k(:,iNuc) = K/nK;
        E1A_ = mI{iNuc}*nK;
        if directAccumulation
          E1A(1:nNucStates(iNuc),iNuc) = E1A_(:);
        else
          E1A(:,iNuc) = E1A_(idxn(:,iNuc)).';
        end
      end
    else
      E1A = 0;
    end

    % second-order
    E2D = 0;
    if secondOrder
      if highSpin
        x =  D1sq*(4*S*(S+1)-3*(8*mS^2-8*mS+3))...
          - D2sq/4*(2*S*(S+1)-3*(2*mS^2-2*mS+1));
        E2D = -x./(2*E0);
      else
        E2DA = 0;
      end
      if (nNuclei>0)
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
          if directAccumulation
            E2A(1:nNucStates(n),n) = E2A_(:);
          else
            E2A(:,n) = E2A_(idxn(:,n)).';
          end
          if highSpin
            DA = Du.'*Ak - uDu*kAu;
            y = DA*(3-6*mS)*mI{n};
            E2DA_ = -y./E0;
            if directAccumulation
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
    
    if directAccumulation
      B0 = (E0-E1D-E2D)*preOri*1e3; % mT
      % compute B shifts
      Bshifts = (-(E1A+E2A+E2DA))*preOri*1e3; % mT
      % directly accumulate into spectrum
      spec = spec + Intensity(imS,iOri)*Exp.AccumWeights(iOri)*...
        multinucstick(B0,nNucStates,Bshifts,Baxis(1),dB,Exp.nPoints);
    else
      if secondOrder
        Bfinal{imS}(iOri,:) = (E0-E1D-E2D-sum(E1A+E2A+E2DA,2))*preOri;
      else
        Bfinal{imS}(iOri,:) = (E0-E1D-sum(E1A,2))*preOri;
      end
    end
    
  end
  
end

if directAccumulation
  B = [];
  Int = [];
  Wid = [];
  Transitions = [];
  spec = spec/dB/prod(nNucStates);
  spec = spec*(2*pi); % powder chi integral
else
  % Positions
  %-------------------------------------------------------------------
  B = [];
  for iTrans = 1:nTransitions
    B = [B Bfinal{iTrans}];
  end
  B = B.'*1e3; % T -> mT
  
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
    lw = sqrt(lw2)*1e6*planck./geff/bmagn*1e3;
    Wid = repmat(lw,nNucTrans*2*S,1);
    
  elseif any(Sys.gStrain) || any(Sys.AStrain)
    
    if any(Sys.gStrain)
      gStrainMatrix = diag(Sys.gStrain./Sys.g)*Exp.mwFreq*1e3; % -> MHz
      if isfield(Sys,'gpa') && any(Sys.gpa(1,:))
        Rp = erot(Sys.gpa(1,:));
        gStrainMatrix = Rp*gStrainMatrix*Rp.';
      end
    else
      gStrainMatrix = zeros(3);
    end
    
    if any(Sys.AStrain) && (Sys.nNuclei>0)
      AStrainMatrix = diag(Sys.AStrain);
      if isfield(Sys,'Apa')
        Rp = erot(Sys.Apa(1,:));
        AStrainMatrix = Rp*AStrainMatrix*Rp.';
      end
      corr = Sys.gAStrainCorr;
      mI1 = -Sys.I(1):+Sys.I(1);
      for idx = 1:numel(mI1)
        StrainMatrix = gStrainMatrix + corr*(mI1(idx))*AStrainMatrix;
        for iOri = 1:nOrientations
          lw2(idx,iOri) = vecs(:,iOri).'*StrainMatrix.^2*vecs(:,iOri);
        end
      end
      Wid = sqrt(lw2)*planck*1e6./repmat(geff,numel(mI1),1)/bmagn*1e3; % MHz -> mT
      idx = repmat(1:numel(mI1),2*S*nNucTrans/numel(mI1),1);
      Wid = Wid(idx(:),:);
    else
      StrainMatrix = gStrainMatrix;
      for iOri = 1:nOrientations
        lw2(1,iOri) = vecs(:,iOri).'*StrainMatrix.^2*vecs(:,iOri);
      end
      Wid = sqrt(lw2)*planck*1e6./geff/bmagn*1e3; % MHz -> mT
      Wid = repmat(Wid,2*S*nNucTrans,1);
    end

  elseif any(Sys.DStrain)
    x = vecs(1,:);
    y = vecs(2,:);
    z = vecs(3,:);
    mS = S:-1:-S;
    mSS = mS.^2-S*(S+1)/3;
    for k = 1:numel(mS)
      dBdD_(k,:) = (3*z.^2-1)/2*mSS(k)*planck./geff/bmagn*1e9;
      dBdE_(k,:) = 3*(x.^2-y.^2)/2*mSS(k)*planck./geff/bmagn*1e9;
    end
    for k = 1:numel(mS)-1
      lwD(k,:) = (dBdD_(k+1,:)-dBdD_(k,:))*Sys.DStrain(1);
      lwE(k,:) = (dBdE_(k+1,:)-dBdE_(k,:))*Sys.DStrain(2);
    end
    lw = sqrt(lwD.^2+lwE.^2);
    Wid = repmat(lw,nNucTrans,1);
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
if (nSites>1) && ~isfield(Opt,'peppercall')
  siz = [nTransitions*nSites, numel(B)/nTransitions/nSites];
  B = reshape(B,siz);
  if ~isempty(Int), Int = reshape(Int,siz); end
  if ~isempty(Wid), Wid = reshape(Wid,siz); end
end

% Arrange output
%---------------------------------------------------------------
Output = {B,Int,Wid,Transitions,spec};
varargout = Output(1:max(nargout,1));

return

%==================================================================

function Combs = allcombinations(varargin)

if (nargin==0), Combs = []; return; end

Combs = varargin{1}(:);
nCombs = numel(Combs);

for iArg = 2:nargin
  New = varargin{iArg}(:);
  nNew = numel(New);
  [idxNew,idxCombs] = find(ones(nNew,nCombs));
  Combs = [Combs(idxCombs(:),:), New(idxNew(:))];
  nCombs = nCombs*nNew;
end

return
