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
DefaultExp.mwFreq = NaN;
DefaultExp.Range = NaN;
DefaultExp.CenterSweep = NaN;
DefaultExp.Temperature = NaN;
DefaultExp.Mode = '';
DefaultExp.Polarization = '';

DefaultExp.CrystalOrientation = [];
DefaultExp.CrystalSymmetry = '';
DefaultExp.MolFrame = [];

Exp = adddefaults(Exp,DefaultExp);

if isnan(Exp.mwFreq), error('Experiment.mwFreq is missing!'); end

if ~isnan(Exp.CenterSweep)
  if ~isnan(Exp.Range)
    %logmsg(0,'Using Experiment.CenterSweep and ignoring Experiment.Range.');
  end
  Exp.Range = Exp.CenterSweep(1) + [-1 1]*Exp.CenterSweep(2)/2;
  Exp.Range = max(Exp.Range,0);
end

if isfield(Exp,'SearchRange'), Exp.Range = Exp.SearchRange; end

if isnan(Exp.Range), error('Experiment.Range/Exp.CenterSweep is missing!'); end
if (diff(Exp.Range)<=0) | ~isfinite(Exp.Range) | ~isreal(Exp.Range) | any(Exp.Range<0)
  error('Exp.Range is not valid!');
end

% Determine excitation mode
p_excitationgeometry;

% Temperature, non-equilibrium populations
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

% Process crystal orientations, crystal symmetry, and frame transforms
% This sets Orientations, nOrientations, nSites and AverageOverChi
p_crystalorientations;

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
c = bmagn/2 * sqrt(S*(S+1)-mS.*(mS-1));
c = c/planck/1e9;
c2 = c.^2;

% Loop over all orientations
for iOri = nOrientations:-1:1
  R = erot(Orientations(iOri,:));
  n0 = R.'*nB0;  % transform to molecular frame representation
  vecs(:,iOri) = n0;
  
  geff(iOri) = norm(g.'*n0);
  u = g.'*n0/geff(iOri); % molecular frame representation
  
  % frequency to field conversion factor
  preOri = 1e6*planck/(geff(iOri)*bmagn);
  
  % Compute intensities
  %----------------------------------------------------------------

  % Compute quantum-mechanical transition rate
  if (linearpolarizedMode)
    if (AverageOverChi)
      TransitionRate(:,iOri) = c2/2*(1-xi1^2)*(trgg-norm(g*u)^2);
    else
      nB1_ = R.'*nB1; % transform to molecular frame representation
      TransitionRate(:,iOri) = c2*(norm(nB1_.'*g)^2-(u.'*g.'*nB1_)^2);
    end
  elseif (unpolarizedMode)
    if (AverageOverChi)
      TransitionRate(:,iOri) = c2/4*(1+xik^2)*(trgg-norm(g*u)^2);
    else
      nk_ = R.'*nk; % transform to molecular frame representation
      TransitionRate(:,iOri) = c2/2*(trgg-norm(g*u)^2-norm(cross(g.'*nk_,u))^2);
    end
  elseif (circpolarizedMode)
    %j = [u(2);-u(1);0]; j = j/norm(j);
    %i = cross(j,u);
    %gixgj = cross(g*i,g*j);
    if (AverageOverChi)
      %TransitionRate(:,iOri) = c2/4*((1+xik^2)*(trgg-norm(g*u)^2)-circpolarizedMode*4*xik^2*n0.'*gixgj);
      TransitionRate(:,iOri) = c2/4*((1+xik^2)*(trgg-norm(g*u)^2)-circpolarizedMode*4*xik^2*det(g)/norm(g.'*n0));
    else
      nk_ = R.'*nk; % transform to molecular frame representation
      %TransitionRate(:,iOri) = c2/2*(trgg-norm(g*u)^2-norm(cross(g.'*nk_,u))^2-circpolarizedMode*2*nk_.'*gixgj);
      TransitionRate(:,iOri) = c2/2*(trgg-norm(g*u)^2-norm(cross(g.'*nk_,u))^2-circpolarizedMode*2*det(g)*xik/norm(g.'*n0));
    end
  end
  if abs(TransitionRate)<1e-10, TransitionRate = 0; end

  % Compute Aasa-Vänngård 1/g factor (frequency-to-field conversion factor)
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
    lw2 = sum(Sys.HStrain.^2*vecs.^2,1); % MHz^2
    lw = sqrt(lw2)*1e6*planck./geff/bmagn*1e3; % mT
    Wid = repmat(lw,nNucTrans*2*S,1);
  else
    Wid = 0;
  end
  
  if any(Sys.gStrain) || any(Sys.AStrain)
    
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
      Wid_gA = sqrt(lw2)*planck*1e6./repmat(geff,numel(mI1),1)/bmagn*1e3; % MHz -> mT
      idx = repmat(1:numel(mI1),2*S*nNucTrans/numel(mI1),1);
      Wid = sqrt(Wid_gA(idx(:),:).^2 + Wid.^2);
    else
      StrainMatrix = gStrainMatrix;
      for iOri = 1:nOrientations
        lw2(1,iOri) = vecs(:,iOri).'*StrainMatrix.^2*vecs(:,iOri);
      end
      Wid_gA = sqrt(lw2)*planck*1e6./geff/bmagn*1e3; % MHz -> mT
      Wid = sqrt(repmat(Wid_gA.^2,2*S*nNucTrans,1)+Wid.^2);
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
    Wid_D = repmat(lw,nNucTrans,1);
    Wid = Sqrt(Wid_D.^2 + Wid.^2);
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
