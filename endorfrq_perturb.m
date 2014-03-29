function varargout = endorfrq_perturb(Sys,Exp,Opt)

% Compute nuclear resonance frequencies based on formulas from
%   Iwasaki, J.Magn.Reson. 16, 417-423 (1974)
% See also
%   Weil, J. Magn. Reson. 113-116 (1975)
%   Hagston and Holmes, J. Phys B: Atom. Molec. Phys. 13, 3505-3519 (1980)

% Assert correct Matlab version
error(chkmlver);

% Check number of input arguments.
switch (nargin)
  case 2, Opt = struct;
  case 3,
  otherwise
    error('Incorrect number of inputs!');
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
if Sys.nElectrons~=1
  err = 'Perturbation theory available only for systems with 1 electron.';
end
if Sys.S~=1/2
%  err = 'Perturbation theory available only for systems with S=1/2.';
end
if any(Sys.AStrain)
  err = ('A strain not supported with perturbation theory. Use matrix diagonalization.');
end
if any(Sys.DStrain)
  err = ('D strain not supported with perturbation theory. Use matrix diagonalization.');
end
error(err);

S = Sys.S;
highSpin = (S>1/2);
SS1 = S*(S+1);
mS = -S:1:S;
mS2 = mS.^2;

if Sys.fullg
  g = Sys.g;
else
  Rg = erot(Sys.gpa);
  g = Rg*diag(Sys.g)*Rg';
end

if highSpin
  if isfield(Sys,'D')
    RD = erot(Sys.Dpa);
    D = diag(Sys.D);
    D = RD*D*RD';
  end
end

I = Sys.I;
nNuclei = Sys.nNuclei;

nuI = -Sys.gn*nmagn*Exp.Field/1e3/planck/1e6;

for iNuc = nNuclei:-1:1
  mI{iNuc} = -I(iNuc):I(iNuc);

  if Sys.fullA
    A{iNuc} = Sys.A((iNuc-1)*3+(1:3),:)*Sys.Ascale(iNuc);
  else
    RA = erot(Sys.Apa(iNuc,:));
    A_ = diag(Sys.A(iNuc,:)*Sys.Ascale(iNuc));
    A{iNuc} = RA*A_*RA.';
  end

  if Sys.fullQ
    P{iNuc} = Sys.Q((iNuc-1)*3+(1:3),:)*Sys.Qscale(iNuc);
  else
    RP = erot(Sys.Qpa(iNuc,:));
    Q_ = diag(Sys.Q(iNuc,:)*Sys.Qscale(iNuc));
    P{iNuc} = RP*Q_*RP.';
  end
  P{iNuc} = P{iNuc} - eye(3)*trace(P{iNuc})/3; % make traceless (for Iwasaki equations)

end

% Experiment
%------------------------------------------------------
err = '';
%if ~isfield(Exp,'mwFreq'), err = 'Exp.mwFreq is missing.'; end
if isfield(Exp,'Detection')
  error('Exp.Detection is obsolete. Use Exp.Mode instead.');
end

if isfield(Exp,'Mode')
  if strcmp(Exp.Mode,'parallel')
    err = 'Parallel mode ENDOR not supported with perturbation theory. Use matrix diagonalization.';
  end
end
if isfield(Exp,'Temperature')
  if (Exp.Temperature<inf)
    err = 'ENDOR temperature effects are not supported with perturbation theory. Use matrix diagonalization.';
  end
end
error(err);

if ~isfield(Exp,'ExciteWidth'), Exp.ExciteWidth = inf; end
OrientationSelection = ~isinf(Exp.ExciteWidth);
lwExcite2 = Exp.ExciteWidth^2;

if ~isfield(Exp,'Orientations')
  logmsg(0,'Exp.Orientations is missing, assuming [0;0]');
  Exp.Orientations = [0;0];
end
Orientations = Exp.Orientations;

% Options
%----------------------------------------------------------------
if isfield(Opt,'Transitions')
  if ~isempty(Opt.Transitions)
    error('ENDOR perturbation theory does not support Opt.Transitions.');
  end
end

if ~isfield(Opt,'PerturbOrder'), Opt.PerturbOrder = 2; end
secondOrder = (Opt.PerturbOrder==2);
if secondOrder
  logmsg(1,'  2nd order perturbation theory');
else
  logmsg(1,'  1st order perturbation theory');
end
if (Opt.PerturbOrder~=1)&&(Opt.PerturbOrder~=2)
  error('Opt.PerturbOrder must be either 1 or 2.');
end

if ~isfield(Opt,'Enhancement'), Opt.Enhancement = 0; end
if (Opt.Enhancement)
  logmsg(1,'  skipping hyperfine enhancement.');
end
logmsg(1,'  hyperfine enhancement: off.');

if ~isfield(Opt,'Nuclei'), Opt.Nuclei = []; end
if isempty(Opt.Nuclei), Opt.Nuclei = 1:nNuclei; end
if any(Opt.Nuclei<1) || any(Opt.Nuclei>nNuclei)
  error('Opt.Nuclei is out of range.');
end


% Orientations
%-----------------------------------------------------------------
[n1,n2] = size(Orientations);
if ((n2==2)||(n2==3)) && (n1~=2) && (n1~=3)
  Orientations = Orientations.';
end
[nAngles,nOrientations] = size(Orientations);
switch nAngles
  case 2,
    AverageOverChi = true;
    Orientations(3,end) = 0; % Entire chi column is set to 0.
  case 3,
    AverageOverChi = false;
  otherwise
    error('Orientations array has %d rows instead of 2 or 3.',nAngles);
end


% Add symmetry-related sites
if ~isfield(Exp,'CrystalSymmetry'), Exp.CrystalSymmetry = ''; end
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
  %nSites = 1;
end


%----------------------------------------------------------
hasQuadrupole = (Sys.I>1/2);

for iNuc = 1:nNuclei
  A_ = A{iNuc};
  detA(iNuc) = det(A_);
  invA{iNuc} = inv(A_);
  trAA(iNuc) = trace(A_'*A_);
end
II1 = I.*(I+1);

if highSpin
  trDD = trace(D^2);
end

states = allcombinations(mI{:});

% Initialize parameters for orientation selectivity determination.
% Selectivity = (maxEPRfreq-minEPRfreq)/minExWidth
if (OrientationSelection)
  maxEPRfreq = -inf;
  minEPRfreq = inf;
  minExWidth = inf;
end

% Loop over all orientations
for iOri = nOrientations:-1:1
  [h1x,h1y,h] = erot(Orientations(:,iOri));
  h = h.';
  
  % zero-order energy: electron Zeeman term only
  geff = norm(g*h);
  E0 = geff*bmagn*Exp.Field*1e-3/planck/1e6;
  u = g*h/geff;
  
  % Compute energies --------------------------------------
  % first- and second-oder zero-field terms
  if highSpin
    Du = D*u;
    uDu = u.'*Du;
    E1D_ = -uDu/2*(SS1-3*mS2);
    if secondOrder
      uDDu = Du.'*Du;
      D1sq = uDDu - uDu^2;
      D2sq = 2*trDD + uDu^2 - 4*uDDu;
      E2D_ = D1sq*mS.*(4*SS1-8*mS2-1) ...
        -D2sq/4*mS.*(2*SS1-2*mS2-1);
      E2D_ = -E2D_/(2*E0);
    else
      E2D_ = zeros(1,numel(mS));
    end
  else
    E1D_ = [0 0];
    E2D_ = [0 0];
  end
  ElectronEnergies = mS*E0 + E1D_ + E2D_;
  
  for iNuc = 1:nNuclei
    mI_ = mI{iNuc};
    Ah = A{iNuc}*g/geff*h; % hyperfine field vector
    for imS = 1:2*S+1
      mS_ = mS(imS);
      
      % first-order hyperfine/nuclear Zeeman term
      HyperfineField = mS_*Ah;
      LarmorField = nuI(iNuc)*h;
      Kh = HyperfineField + LarmorField; % total field vector at nucleus
      normKh = norm(Kh);
      k = Kh/normKh;
      E1A_ = mI_*normKh;
      strongCouplingCase = norm(HyperfineField) > norm(LarmorField);
      if strongCouplingCase && (mS_<0), E1A_ = -E1A_; end

      % first-order nuclear quadrupole term
      if hasQuadrupole(iNuc)
        kPk = k.'*P{iNuc}*k;
        E1Q_ = -kPk/2*(II1(iNuc)-3*mI_.^2);
      else
        E1Q_ = 0;
      end

      if secondOrder
        
        % second-order hyperfine term
        Ak = A{iNuc}.'*k;
        kAu = Ak.'*u;
        kAAk = Ak.'*Ak;
        A1sq = kAAk - kAu^2;
        A2 = detA(iNuc)*(u.'*invA{iNuc}*k);
        A3 = trAA(iNuc) - norm(A{iNuc}*u)^2 - kAAk + kAu^2;
        %E2A_ =  A1sq*mS_*mI_.^2 ...
        %     - A2*(SS1-mS_^2)*mI_*sign(mS_) ... % why sign?? seems wrong
        %     + A3/2*mS_*(II1(n)-mI_.^2);
        E2A_ =  A1sq*mS_*mI_.^2 ...
             - A2*(SS1-mS_^2)*mI_ ...
             + A3/2*mS_*(II1(iNuc)-mI_.^2);
        E2A_ = E2A_/(2*E0);
        
        % second-order zerofield-hyperfine crossterm
        if highSpin
          DA = Du.'*Ak - uDu*kAu;
          E2DA = DA*(SS1-3*mS_^2)*mI_;
          E2DA = -E2DA/E0;
        else
          E2DA = 0;
        end
        
        % second-order nuclear quadrupole term
        if hasQuadrupole(iNuc)
          kPPk = k.'*P{iNuc}^2*k;
          P1sq = kPPk - kPk^2;
          P2sq = 2*trace(P{iNuc})^2 + kPk^2 - 4*kPPk;
          E2Q_ = P1sq*mI_.*(4*II1(iNuc)-8*mI_.^2-1)...
              -P2sq/4*mI_.*(2*II1(iNuc)-2*mI_.^2-1);
          E2Q_ = -E2Q_/(2*normKh);
        else
          E2Q_ = 0;
        end
        
        NuclearEnergyShifts{iNuc}(imS,:) = E1Q_+E1A_+E2A_+E2DA+E2Q_;

      else
        
        NuclearEnergyShifts{iNuc}(imS,:) = E1Q_+E1A_;
      
      end % if secondOrder

    end % for imS = 1:2*S+1
  end % for n = 1:nNuclei
  % --------------------------------------------------------
  
  % ENDOR frequencies --------------------------------------
  for iNuc = 1:nNuclei
    dE_ = abs(diff(NuclearEnergyShifts{iNuc},1,2));
    EndorFrequencies{iNuc}(iOri,:) = dE_(:).';
  end
  % --------------------------------------------------------

  
  % Orientation selection ----------------------------------
  if (OrientationSelection)
    
    % compute EPR transition frequencies
    for imS = 1:2*S+1
      for iNuc = 1:nNuclei
        EnergyShifts(:,iNuc) = NuclearEnergyShifts{iNuc}(imS,states(:,iNuc));
      end
      TotalEnergy(:,imS) = ElectronEnergies(imS) + sum(EnergyShifts,2);
    end
    EPRfreq = diff(TotalEnergy,1,2);
    
    % compute excitation factors for EPR transitions
    lw2 = sum((Sys.HStrain.*h.').^2) + lwExcite2;
    offsetFreq = EPRfreq - Exp.mwFreq*1e3;
    exciteFactor = exp(-2/lw2*offsetFreq.^2);

    % collect statistics
    if (min(EPRfreq(:))<minEPRfreq), minEPRfreq = min(EPRfreq(:)); end
    if (max(EPRfreq(:))>maxEPRfreq), maxEPRfreq = max(EPRfreq(:)); end
    if (lw2<minExWidth^2), minExWidth = sqrt(lw2); end

    % for each ENDOR transition, combine excitation factors of
    % EPR transitions that share a level with the ENDOR transition
    for iNuc = 1:nNuclei
      for imI = 1:numel(mI{iNuc})-1
        idx = (states(:,iNuc)==imI) | (states(:,iNuc)==imI+1);
        totalWeight{iNuc}(iOri,imI) = sum(exciteFactor(idx));
      end
    end
    % totalWeight{iNuc} size: #orientations x (2*I+1)
  end
  %---------------------------------------------------------
  
end % for iOri = ...

% Compute selectivity.
if (OrientationSelection)
  Selectivity = (maxEPRfreq-minEPRfreq)/minExWidth;
  logmsg(2,'  limited excitation width, selectivity %g',Selectivity);
else
  Selectivity = 0;
  logmsg(2,'  infinite excitation width');
end
Info.Selectivity = Selectivity;

% Frequencies
%-------------------------------------------------------------------
EndorFreqs = [];
for iNuc = Opt.Nuclei
  EndorFreqs = [EndorFreqs; EndorFrequencies{iNuc}.'];
end

% Intensities
%-------------------------------------------------------------------
EndorIntensities = [];
nNucStates = 2*I+1;
for iNuc = Opt.Nuclei
  matrixelement = sqrt((1:2*I(iNuc)).*(2*I(iNuc):-1:1))/2;
  pre = Sys.gn(iNuc)*nmagn/planck/1e6/1e3; % MHz/mT
  pre = (pre*matrixelement).^2;
  if AverageOverChi, end
  newIntensities = pre(:)/nNucStates(iNuc);
  newIntensities = repmat(newIntensities,1,nOrientations);
  % array size now: (2*I(iNuc)) x nOrientations
  if (OrientationSelection)
    newIntensities = newIntensities.*totalWeight{iNuc}.';
  end
  idx = 1:2*I(iNuc); idx = repmat(idx,2*S+1,1); idx = idx(:);
  newIntensities = newIntensities(idx,:);
  % transitions(=rows) ordered by 1) mI of lower level, 2) mS
  % e.g. (-I,+1/2), (-I,-1/2), (-I+1,+1/2), (-I+1,-1/2), etc
  % (identical to ordering of EndorFreqs).
  EndorIntensities = [EndorIntensities; newIntensities];
  
end
EndorIntensities = EndorIntensities*prod(2*Sys.I+1);

% Arrange output
%---------------------------------------------------------------
Transitions = [];
Output = {EndorFreqs,EndorIntensities,Transitions,Info};
varargout = Output(1:max(nargout,1));

function Combs = allcombinations(varargin)

if (nargin==0), Combs = []; return; end

idxonly = 1;
if idxonly
  Combs = (1:numel(varargin{1})).';
else
  Combs = varargin{1}(:);
end
nCombs = numel(Combs);

for iArg = 2:nargin
  New = varargin{iArg}(:);
  nNew = numel(New);
  [idxNew,idxCombs] = find(ones(nNew,nCombs));
  if idxonly
    Combs = [Combs(idxCombs(:),:), idxNew(:)];
  else
    Combs = [Combs(idxCombs(:),:), New(idxNew(:))];
  end
  nCombs = nCombs*nNew;
end

return
