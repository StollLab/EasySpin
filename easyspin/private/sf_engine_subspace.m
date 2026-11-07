% sf_engine_subspace  Subspace engine of saffron
%
%   [y,info] = sf_engine_subspace(Sys,Exp,Opt,Orientations)
%
% Simulates a pulse EPR signal as a sum over electron coherence transfer
% pathways, in the nuclear subspaces of the electron spin manifolds. This
% generalizes the formalism of Mims to arbitrary sequences:
%   W. B. Mims, Phys. Rev. B 5, 2409 (1972), doi:10.1103/PhysRevB.5.2409
%
% Input:
%   Sys           validated spin system (single isotopologue)
%   Exp           experiment, with Exp.OriWeights
%   Opt           options, with the parsed sequence in Opt.fastSeq and the
%                 canonical Exp in Opt.fastExp
%   Orientations  orientations (phi,theta,chi), one per row
%
% Output:
%   y     signal, with one array dimension per indirect dimension
%         (Exp.Dim1, Exp.Dim2, ...), normalized and phased (Exp.DetPhase)
%   info  structure with decay rates (for processing) and engine details
%
% Physics (see the saffron documentation and the redesign notes):
% - Each excited EPR transition (b -> a) is a two-level electron system.
%   The nuclear sub-Hamiltonians of the two manifolds use the expectation
%   value <S> of the electron spin (high-field approximation). Their
%   eigenvectors Va and Vb give the overlap matrix M = Va'*Vb.
% - Ideal microwave pulses give pathway prefactors (sf_parsesequence);
%   rectangular pulses are exact 2n x 2n propagators at each resonance
%   offset. Free evolution is a diagonal phase factor.
% - Linearly incremented delays: the signal of each pathway is a sum of
%   exponentials, binned in the frequency domain (sf_peaks) and Fourier
%   transformed.
% - RF pulses (ENDOR, rectangular only): population transfer between the
%   levels of each NMR transition, with the two-level response (Rabi
%   formula) at the detuning from the RF carrier, averaged over the RF phase.
% - Relaxation: exp(-t/T2) for intervals with electron coherence and
%   exp(-t/T1) for intervals with electron populations.

function [y,info] = sf_engine_subspace(Sys,Exp,Opt,Orientations)

seq = Opt.fastSeq;
nOrientations = size(Orientations,1);
info = struct;

% Relaxation times (0 or not given: no relaxation)
T1 = 0; T2 = 0;
if isfield(Sys,'T1') && ~isempty(Sys.T1)
  if numel(Sys.T1)~=1
    error('The subspace engine does not support transition-selective relaxation times (Sys.T1). Use Opt.Method = ''fullspace''.');
  end
  T1 = Sys.T1;
end
if isfield(Sys,'T2') && ~isempty(Sys.T2)
  if numel(Sys.T2)~=1
    error('The subspace engine does not support transition-selective relaxation times (Sys.T2). Use Opt.Method = ''fullspace''.');
  end
  T2 = Sys.T2;
end
if T1<0 || T2<0 || ~isreal(T1) || ~isreal(T2)
  error('T1 and T2 in Sys.T1 and Sys.T2 must be positive, in microseconds.');
end
if T1==0, T1 = Inf; end
if T2==0, T2 = Inf; end
if isfield(Opt,'Relaxation') && ~isempty(Opt.Relaxation) && ~any(Opt.Relaxation)
  T1 = Inf; T2 = Inf;
end
relaxation = ~isinf(T1) || ~isinf(T2);

isENDOR = any(seq.kind==2);

% Nuclei included in the nuclear subspaces
%-------------------------------------------------------------------------------
if isfield(Opt,'Nuclei') && ~isempty(Opt.Nuclei)
  if any(Opt.Nuclei>Sys.nNuclei) || any(Opt.Nuclei<1)
    error('Index in Opt.Nuclei out of range.');
  end
  shfNuclei = Opt.Nuclei;
else
  shfNuclei = 1:Sys.nNuclei;
  % With ideal pulses and an excitation window (Opt.IdealPulses), nuclei
  % with hyperfine couplings larger than the window are moved to the core
  % system: they only shift the EPR transitions. With finite pulses, all
  % nuclei are included and the pulses select the transitions.
  if ~isempty(shfNuclei) && ~isempty(seq.ExciteWidth) && ~seq.finitePulses
    strong = false(1,Sys.nNuclei);
    for iNuc = 1:Sys.nNuclei
      if Sys.fullA
        A_ = Sys.A((iNuc-1)*3+(1:3),:);
      else
        A_ = Sys.A(iNuc,:);
      end
      strong(iNuc) = max(abs(A_(:)))>seq.ExciteWidth;
    end
    shfNuclei(strong) = [];
  end
end
twoElectronManifolds = Sys.nElectrons==1 && Sys.S==1/2 && numel(shfNuclei)==Sys.nNuclei;

% Product rule: exact for ideal pulses and S = 1/2, but it needs a Fourier
% transform per entry, pathway and nucleus. Benchmarks (2026-11-03) found
% it slower than the full nuclear space in most cases, and it adds binning
% errors, so it is used only on request (Opt.ProductRule = true).
if ~isfield(Opt,'ProductRule') || isempty(Opt.ProductRule) || strcmp(Opt.ProductRule,'auto')
  productRule = false;
else
  productRule = logical(Opt.ProductRule);
  if productRule && (seq.finitePulses || isENDOR)
    error('Opt.ProductRule requires ideal pulses and no RF pulses.');
  end
end
if numel(shfNuclei)<=1, productRule = false; end
info.ProductRule = productRule;

% Nuclear sub-Hamiltonians
%-------------------------------------------------------------------------------
NucHams = nuclearHamiltonians(Sys,shfNuclei,productRule);
nSubSpaces = numel(NucHams);

% Electronic core system (not needed for S = 1/2 with all nuclei)
if ~twoElectronManifolds
  coreSys = nucspinrmv(Sys,shfNuclei);
  coreSys.processed = false;
  coreSys.lw = 0;
  if isfield(coreSys,'A_'), coreSys = rmfield(coreSys,'A'); end
  [H0,mux,muy,muz] = ham(coreSys);
  Sx = sop(coreSys,[1,1]);
  Sy = sop(coreSys,[1,2]);
  Sz = sop(coreSys,[1,3]);
else
  g = Sys.g;
  if ~Sys.fullg
    Rg = erot(Sys.gFrame).';
    g = Rg*diag(g)*Rg.';
  end
end

OrientationSelection = isfield(Exp,'mwFreq') && ~isempty(seq.ExciteWidth);
if ~isfield(Opt,'OriThreshold'), Opt.OriThreshold = 0.005; end
if ~isfield(Sys,'lwAnisoG'), Sys.lwAnisoG = [0 0 0]; end
lwAnisoG = Sys.lwAnisoG(:).'.*[1 1 1];

% Field-reversed orientations (partial-sphere powder grids with nuclear
% quadrupole interaction, see the symmetry probe below)
partialSphereGrid = ~Opt.crystalSample && isfield(Opt,'nOctants') && Opt.nOctants~=8;
hasQ = nSubSpaces>0 && any(arrayfun(@(N)any(N.Hnq(:)~=0),NucHams));
reversedField = "none";
if partialSphereGrid && hasQ
  reversedField = "compute";
  if ~seq.finitePulses && ~isENDOR
    reversedField = symmetryProbe(seq);
  end
end
info.reversedField = reversedField;

% Entries: one per excited EPR transition and orientation (and reversed
% field), with nuclear energies and overlap matrices of all subspaces
%-------------------------------------------------------------------------------
E = struct;
E.w = zeros(1,0);            % weight (orientation weight x selection weight)
E.dE0 = zeros(1,0);          % EPR resonance offset from mwFreq (MHz)
E.lw = zeros(1,0);           % inhomogeneous FWHM (MHz)
E.ori = zeros(1,0);
E.rev = false(1,0);
for s = 1:nSubSpaces
  E.Ea{s} = []; E.Eb{s} = []; E.M{s} = [];
  if isENDOR, E.Va{s} = []; E.Vb{s} = []; end
end
if twoElectronManifolds
  zL = zeros(3,nOrientations); xL = zL; yL = zL;
  for iOri = 1:nOrientations
    [xL(:,iOri),yL(:,iOri),zL(:,iOri)] = erot(Orientations(iOri,:),'rows');
  end
  geff = g*zL;
  geffNorm = sqrt(sum(geff.^2,1));
  nuEPR = geffNorm*bmagn*Exp.Field/1e3/planck/1e6;  % MHz
  lwInh = sqrt(sum((lwAnisoG(:).*zL).^2,1));
  if OrientationSelection
    lw2 = lwInh.^2 + seq.ExciteWidth^2;
    oriSel = exp(-4*log(2)*(nuEPR-1000*Exp.mwFreq).^2./lw2)./sqrt(1+lwInh.^2/seq.ExciteWidth^2);
  else
    oriSel = ones(1,nOrientations);
  end
  keepOri = oriSel>=Opt.OriThreshold;
  info.nSkippedOrientations = sum(~keepOri);
  oriList = find(keepOri);
  nE = numel(oriList);
  E.w = Exp.OriWeights(oriList).*oriSel(oriList);
  if isfield(Exp,'mwFreq')
    E.dE0 = nuEPR(oriList) - 1000*Exp.mwFreq;
  else
    E.dE0 = zeros(1,nE);
  end
  E.lw = lwInh(oriList);
  E.ori = oriList;
  E.rev = false(1,nE);
  E.a = 2*ones(1,nE); E.b = ones(1,nE);
  for s = 1:nSubSpaces
    n = size(NucHams(s).Hnzz,1);
    Ea = zeros(n,nE); Eb = zeros(n,nE); M = zeros(n,n,nE);
    if isENDOR, VaAll = M; VbAll = M; end
    for k = 1:nE
      iOri = oriList(k);
      q = geff(:,iOri).'/geffNorm(iOri);
      Hnuc = Exp.Field*(zL(1,iOri)*NucHams(s).Hnzx + zL(2,iOri)*NucHams(s).Hnzy + zL(3,iOri)*NucHams(s).Hnzz) + NucHams(s).Hnq;
      % manifold a: <S> = +q/2, manifold b: <S> = -q/2
      Ha = Hnuc + 0.5*(q(1)*NucHams(s).Hhfx + q(2)*NucHams(s).Hhfy + q(3)*NucHams(s).Hhfz);
      Hb = Hnuc - 0.5*(q(1)*NucHams(s).Hhfx + q(2)*NucHams(s).Hhfy + q(3)*NucHams(s).Hhfz);
      [Va,Ea_] = eig((Ha+Ha')/2,'vector');
      [Vb,Eb_] = eig((Hb+Hb')/2,'vector');
      Ea(:,k) = real(Ea_); Eb(:,k) = real(Eb_);
      M(:,:,k) = Va'*Vb;
      if isENDOR, VaAll(:,:,k) = Va; VbAll(:,:,k) = Vb; end
    end
    E.Ea{s} = Ea; E.Eb{s} = Eb; E.M{s} = M;
    if isENDOR, E.Va{s} = VaAll; E.Vb{s} = VbAll; end
  end
  if isENDOR
    E.geff = geff(:,oriList); E.geffNorm = geffNorm(oriList);
    E.xL = xL(:,oriList); E.yL = yL(:,oriList);
  end
else
  % General electron system: transitions from the core Hamiltonian
  rows = {};
  for iOri = 1:nOrientations
    [~,yLab_M,zLab_M] = erot(Orientations(iOri,:),'rows');
    muzL = zLab_M(1)*mux + zLab_M(2)*muy + zLab_M(3)*muz;
    H = H0 - Exp.Field*muzL;
    [eV,eE] = eig(H);
    eE = real(diag(eE));
    SyLab = yLab_M(1)*Sx + yLab_M(2)*Sy + yLab_M(3)*Sz;
    SyLab = abs(eV'*SyLab*eV);
    maxSyy = max(SyLab(:));
    lwInh2 = sum((lwAnisoG(:).*zLab_M(:)).^2);
    if OrientationSelection
      lw2 = lwInh2 + seq.ExciteWidth^2;
      dE = eE - eE.' - Exp.mwFreq*1e3;
      SyLab = SyLab.*exp(-4*log(2)*dE.^2/lw2)/sqrt(1+lwInh2/seq.ExciteWidth^2);
    end
    SyLab(SyLab<Opt.OriThreshold*maxSyy) = 0;
    if isfield(Opt,'Transitions') && ~isempty(Opt.Transitions)
      nStates = size(H,1);
      if any(Opt.Transitions(:)>nStates)
        error('Opt.Transitions out of range - only values between 1 and %d are valid.',nStates);
      end
      rmvTransition = true(size(H));
      for tr_ = 1:size(Opt.Transitions,1)
        tr = Opt.Transitions(tr_,:);
        rmvTransition(tr(1),tr(2)) = false;
        rmvTransition(tr(2),tr(1)) = false;
      end
      SyLab(rmvTransition) = 0;
    end
    [v,u,wSel] = find(tril(SyLab,-1));
    for iT = 1:numel(u)
      b = u(iT); a = v(iT);
      vecA = eV(:,a); vecB = eV(:,b);
      SA = real([vecA'*Sx*vecA; vecA'*Sy*vecA; vecA'*Sz*vecA]);
      SB = real([vecB'*Sx*vecB; vecB'*Sy*vecB; vecB'*Sz*vecB]);
      rows(end+1,:) = {iOri, Exp.OriWeights(iOri)*wSel(iT), eE(a)-eE(b), sqrt(lwInh2), SA, SB, zLab_M}; %#ok<AGROW>
    end
  end
  nE = size(rows,1);
  info.nSkippedOrientations = nOrientations - numel(unique([rows{:,1}]));
  E.w = [rows{:,2}];
  if isfield(Exp,'mwFreq')
    E.dE0 = [rows{:,3}] - 1000*Exp.mwFreq;
  else
    E.dE0 = zeros(1,nE);
  end
  E.lw = [rows{:,4}];
  E.ori = [rows{:,1}];
  E.rev = false(1,nE);
  for s = 1:nSubSpaces
    n = size(NucHams(s).Hnzz,1);
    Ea = zeros(n,nE); Eb = zeros(n,nE); M = zeros(n,n,nE);
    if isENDOR, VaAll = M; VbAll = M; end
    for k = 1:nE
      zLab_M = rows{k,7};
      Hnuc = Exp.Field*(zLab_M(1)*NucHams(s).Hnzx + zLab_M(2)*NucHams(s).Hnzy + zLab_M(3)*NucHams(s).Hnzz) + NucHams(s).Hnq;
      SA = rows{k,5}; SB = rows{k,6};
      Ha = Hnuc + SA(1)*NucHams(s).Hhfx + SA(2)*NucHams(s).Hhfy + SA(3)*NucHams(s).Hhfz;
      Hb = Hnuc + SB(1)*NucHams(s).Hhfx + SB(2)*NucHams(s).Hhfy + SB(3)*NucHams(s).Hhfz;
      [Va,Ea_] = eig((Ha+Ha')/2,'vector');
      [Vb,Eb_] = eig((Hb+Hb')/2,'vector');
      Ea(:,k) = real(Ea_); Eb(:,k) = real(Eb_);
      M(:,:,k) = Va'*Vb;
      if isENDOR, VaAll(:,:,k) = Va; VbAll(:,:,k) = Vb; end
    end
    E.Ea{s} = Ea; E.Eb{s} = Eb; E.M{s} = M;
    if isENDOR, E.Va{s} = VaAll; E.Vb{s} = VbAll; end
  end
  if isENDOR
    E.zL = [rows{:,7}];
    E.geff = []; % no enhancement model for the general electron system
  end
end
if nSubSpaces==0
  error('The subspace engine needs at least one nucleus with nuclear spin in the nuclear subspace (Opt.Nuclei).');
end

% Field-reversed partners, computed explicitly (M, Va, Vb -> conj)
if reversedField=="compute"
  E.w = [E.w E.w]/2;
  E.dE0 = [E.dE0 E.dE0];
  E.lw = [E.lw E.lw];
  E.ori = [E.ori E.ori];
  E.rev = [false(1,nE) true(1,nE)];
  for s = 1:nSubSpaces
    E.Ea{s} = [E.Ea{s} E.Ea{s}];
    E.Eb{s} = [E.Eb{s} E.Eb{s}];
    E.M{s} = cat(3,E.M{s},conj(E.M{s}));
    if isENDOR
      E.Va{s} = cat(3,E.Va{s},conj(E.Va{s}));
      E.Vb{s} = cat(3,E.Vb{s},conj(E.Vb{s}));
    end
  end
  if isENDOR
    oriFields = intersect(fieldnames(E),{'geff','geffNorm','xL','yL','zL'});
    for iField = 1:numel(oriFields)
      fld = oriFields{iField};
      if ~isempty(E.(fld)), E.(fld) = [E.(fld) E.(fld)]; end
    end
  end
  nE = 2*nE;
end
E.nE = nE;
info.nEntries = nE;

% Offsets (finite pulses): Gaussian distribution of each entry, FWHM from
% Sys.lwAnisoG, on a uniform grid over +-2 widths
if seq.finitePulses
  anyInh = any(E.lw>0);
  if ~isfield(Opt,'nOffsets') || isempty(Opt.nOffsets)
    if anyInh, Opt.nOffsets = 291; else, Opt.nOffsets = 1; end
  end
  if Opt.nOffsets==1 || ~anyInh
    offsetGrid = 0; offsetWeight = 1;
  else
    offsetGrid = linspace(-2,2,Opt.nOffsets);  % in units of the 1/e half width
    offsetWeight = exp(-offsetGrid.^2);
    offsetWeight = offsetWeight/sum(offsetWeight);
  end
else
  offsetGrid = 0; offsetWeight = 1;
end
E.offsetGrid = offsetGrid;
E.offsetWeight = offsetWeight;

% Simulation
%-------------------------------------------------------------------------------
Rel.T1 = T1; Rel.T2 = T2; Rel.on = relaxation;
expandFactor = 4;  % 2^expandFactor buffer expansion for binning
if numel(seq.kind1Dims)==2, expandFactor = 2; end
if isfield(Opt,'Expand') && ~isempty(Opt.Expand), expandFactor = Opt.Expand; end

if isENDOR
  y = endorSweep(E,seq,Sys,Exp,NucHams,nSubSpaces);
  rates = [];
else
  kind3 = find(seq.kind==3);
  if isempty(kind3)
    [y,rates] = pathwaySignal(E,seq,Rel,productRule,expandFactor);
  else
    % Point-by-point dimensions: re-parse the sequence with modified pulses
    % and delays for each point
    nPts = seq.nPoints(kind3);
    yk = cell([nPts 1]);
    for iPt = 1:prod(nPts)
      sub = cell(1,numel(kind3));
      [sub{:}] = ind2sub([nPts 1],iPt);
      Expk = Opt.fastExp;
      pulseIdx = find(cellfun(@isstruct,Expk.Sequence));
      delayIdx = find(~cellfun(@isstruct,Expk.Sequence));
      for j = 1:numel(kind3)
        d = kind3(j);
        m = seq.mods{d};
        for r = 1:size(m,1)
          val = m{r,4}(sub{j});
          if m{r,1}=='p'
            el = Expk.Sequence{pulseIdx(m{r,2})};
            if ~isfield(el,m{r,3}) || isempty(el.(m{r,3}))
              switch m{r,3}
                case {'Phase','Frequency','tp'}, el.(m{r,3}) = 0;
                otherwise, error('Pulse field %s is not set and cannot be incremented.',m{r,3});
              end
            end
            el.(m{r,3}) = el.(m{r,3}) + val;
            Expk.Sequence{pulseIdx(m{r,2})} = el;
          else
            Expk.Sequence{delayIdx(m{r,2})} = Expk.Sequence{delayIdx(m{r,2})} + val;
          end
        end
      end
      % remove kind-3 dimensions
      keepDims = setdiff(1:numel(seq.kind),kind3);
      Expk.nPoints = seq.nPoints(keepDims);
      for d = 1:numel(seq.kind), Expk = rmfield(Expk,sprintf('Dim%d',d)); end
      for j = 1:numel(keepDims)
        Expk.(sprintf('Dim%d',j)) = Opt.fastExp.(sprintf('Dim%d',keepDims(j)));
      end
      if isempty(Expk.nPoints), Expk = rmfield(Expk,'nPoints'); end
      [seqk,reasons] = sf_parsesequence(Expk,Opt);
      if ~isempty(reasons)
        error('The subspace engine cannot simulate point %d: %s.',iPt,strjoin(reasons,', '));
      end
      [yk{iPt},rates] = pathwaySignal(E,seqk,Rel,productRule,expandFactor);
    end
    % Assemble: kind-1 dims first, then kind-3 dims, then permute
    k1 = seq.kind1Dims;
    sz1 = seq.nPoints(k1);
    Y = zeros(prod(sz1),prod(nPts));
    for iPt = 1:prod(nPts)
      Y(:,iPt) = yk{iPt}(:);
    end
    y = reshape(Y,[sz1 nPts 1 1]);
    perm = [k1 kind3];
    [~,inv] = sort(perm);
    if numel(perm)>1
      y = permute(y,[inv numel(perm)+1:max(3,numel(perm)+1)]);
    end
  end
end
if isvector(y) && ~isscalar(y), y = y(:).'; end
info.decayRates = rates;

% Detection phase, field-reversed partners, normalization
y = y*exp(-1i*seq.DetPhase);
switch reversedField
  case "conj"
    y = real(y);
  case "conjtranspose"
    y = (y + y')/2;
end
% Normalization: true echo amplitude Tr(S+*rho)/Tr(S+'*S+) for the initial
% state -Sz (S = 1/2). For other spin systems, transitions are weighted
% by |<u|Sy|v>| and are not on this scale.
y = y/prod(2*Sys.I+1);

end


%===============================================================================
% Nuclear sub-Hamiltonians (MHz, per mT for the nuclear Zeeman terms)
function NucHams = nuclearHamiltonians(Sys,shfNuclei,productRule)

NucHams = struct('Hnzx',{},'Hnzy',{},'Hnzz',{},'Hhfx',{},'Hhfy',{},'Hhfz',{},'Hnq',{});
if isempty(shfNuclei), return; end
I = Sys.I(shfNuclei);
if productRule
  nSub = numel(shfNuclei);
else
  nSub = 1;
  dim = prod(2*I+1);
  z = zeros(dim);
  NucHams(1) = struct('Hnzx',z,'Hnzy',z,'Hnzz',z,'Hhfx',z,'Hhfy',z,'Hhfz',z,'Hnq',z);
end
for ii = 1:numel(shfNuclei)
  iNuc = shfNuclei(ii);
  if productRule
    Ix = sop(I(ii),[1,1]); Iy = sop(I(ii),[1,2]); Iz = sop(I(ii),[1,3]);
    s = ii;
    z = zeros(size(Ix));
    NucHams(s) = struct('Hnzx',z,'Hnzy',z,'Hnzz',z,'Hhfx',z,'Hhfy',z,'Hhfz',z,'Hnq',z);
  else
    Ix = sop(I,[ii,1]); Iy = sop(I,[ii,2]); Iz = sop(I,[ii,3]);
    s = 1;
  end
  Ix = full(Ix); Iy = full(Iy); Iz = full(Iz);
  pre = -Sys.gn(iNuc)*nmagn/1e3/planck/1e6;  % MHz/mT
  NucHams(s).Hnzx = NucHams(s).Hnzx + pre*Ix;
  NucHams(s).Hnzy = NucHams(s).Hnzy + pre*Iy;
  NucHams(s).Hnzz = NucHams(s).Hnzz + pre*Iz;
  if Sys.fullA
    A = Sys.A((iNuc-1)*3+(1:3),:);
  else
    if isfield(Sys,'AFrame')
      R_A2M = erot(Sys.AFrame(iNuc,:)).';
    else
      R_A2M = eye(3);
    end
    A = R_A2M*diag(Sys.A(iNuc,:))*R_A2M.';
    A = (A+A.')/2;
  end
  NucHams(s).Hhfx = NucHams(s).Hhfx + A(1,1)*Ix + A(1,2)*Iy + A(1,3)*Iz;
  NucHams(s).Hhfy = NucHams(s).Hhfy + A(2,1)*Ix + A(2,2)*Iy + A(2,3)*Iz;
  NucHams(s).Hhfz = NucHams(s).Hhfz + A(3,1)*Ix + A(3,2)*Iy + A(3,3)*Iz;
  if Sys.I(iNuc)>=1
    if Sys.fullQ
      Q = Sys.Q((iNuc-1)*3+(1:3),:);
    else
      if isfield(Sys,'QFrame')
        R_Q2M = erot(Sys.QFrame(iNuc,:)).';
      else
        R_Q2M = eye(3);
      end
      Q = R_Q2M*diag(Sys.Q(iNuc,:))*R_Q2M.';
      Q = (Q+Q.')/2;
    end
    NucHams(s).Hnq = NucHams(s).Hnq + ...
      Ix*(Q(1,1)*Ix + Q(1,2)*Iy + Q(1,3)*Iz) + ...
      Iy*(Q(2,1)*Ix + Q(2,2)*Iy + Q(2,3)*Iz) + ...
      Iz*(Q(3,1)*Ix + Q(3,2)*Iy + Q(3,3)*Iz);
  end
end
end


%===============================================================================
% Signal over the linearly incremented dimensions (kind 1), for all
% entries, offsets and pathways
function [y,rates] = pathwaySignal(E,seq,Rel,productRule,expandFactor)

N = seq.N;
nP = seq.nPathways;
incIntervals = find(seq.Inc~=0);
nDims = numel(seq.kind1Dims);
scheme = seq.IncSchemeID;

% Relaxation: fixed decay factor per pathway and decay rates per dimension;
% pathways with the same rates share a buffer
coh = seq.co~=0;
rateInt = zeros(nP,N);
if Rel.on
  rateInt(coh) = 1/Rel.T2;
  rateInt(~coh) = 1/Rel.T1;
end
fixedDecay = exp(-rateInt*seq.t.');
pRates = zeros(nP,nDims);
for d = 1:nDims
  pRates(:,d) = sum(rateInt(:,seq.Inc==d),2);
end
[rates,~,group] = unique(pRates,'rows');
nGroups = size(rates,1);
pref = seq.pathwayprefactor(:).*fixedDecay;

% Few points: evaluate each point exactly (no frequency binning)
if nDims>0 && prod(seq.nPointsInc)<=64 && ~productRule
  if nDims==1, y = zeros(1,seq.nPointsInc); else, y = zeros(seq.nPointsInc); end
  for iPt = 1:prod(seq.nPointsInc)
    [i1,i2] = ind2sub([seq.nPointsInc 1],iPt);
    k = [i1 i2]-1;
    s = seq;
    for i = incIntervals
      s.t(i) = seq.t(i) + k(seq.Inc(i))*seq.dt(seq.Inc(i));
    end
    s.Inc(:) = 0;
    s.kind1Dims = [];
    s.dt = [];
    s.nPointsInc = [];
    s.IncSchemeID = 0;
    s.idxIncL = zeros(nP,0);
    s.idxIncR = zeros(nP,0);
    y(iPt) = pathwaySignal(E,s,Rel,false,expandFactor);
  end
  return
end

% Buffers
if nDims>0
  nPointsF = 2^expandFactor*seq.nPointsInc;
  if nDims==1, siz = [1 nPointsF]; else, siz = nPointsF; end
end

nSub = numel(E.Ea);
if productRule
  % per pathway and subspace, time domain products, per entry
  y = 0;
  exact = prod(seq.nPointsInc)<=64;
  for e = 1:E.nE
    for p = 1:nP
      td_ = 1;
      for s = 1:nSub
        Ed = entrySlice(E,s,e);
        if exact
          % few points: exact evaluation of each point
          if nDims==1, td1 = zeros(1,seq.nPointsInc); else, td1 = zeros(seq.nPointsInc); end
          for iPt = 1:numel(td1)
            [i1,i2] = ind2sub([seq.nPointsInc 1],iPt);
            k = [i1 i2]-1;
            sk = seq;
            for i = incIntervals
              sk.t(i) = seq.t(i) + k(seq.Inc(i))*seq.dt(seq.Inc(i));
            end
            sk.Inc(:) = 0;
            [~,D] = evaluatePathways(Ed,sk,p,[]);
            td1(iPt) = trace(D);
          end
          if s==1, td1 = td1*pref(p)*E.w(e); end
          td_ = td_.*td1;
          continue
        end
        [G,D,BL,BR] = evaluatePathways(Ed,seq,p,[]);
        if s==1, G = G*pref(p)*E.w(e); end
        bRe = zeros(siz); bIm = zeros(siz);
        sf_peaks(scheme,bRe,bIm,seq.dt,seq.idxIncL(p,:),seq.idxIncR(p,:),...
          Ed.Ea,Ed.Eb,G,D,BL{:},BR{:});
        buff = complex(bRe,bIm);
        if nDims==1
          td1 = ifft(buff)*numel(buff);
          td1 = td1(1:seq.nPointsInc);
        else
          td1 = ifft2(buff)*numel(buff);
          td1 = td1(1:seq.nPointsInc(1),1:seq.nPointsInc(2));
        end
        td_ = td_.*td1;
      end
      td_ = applyDecay(td_,rates(group(p),:),seq);
      y = y + td_;
    end
  end
  return
end

% Without product rule: one subspace
if nDims>0
  bRe = cell(1,nGroups); bIm = cell(1,nGroups);
  for iG = 1:nGroups
    bRe{iG} = zeros(siz); bIm{iG} = zeros(siz);
  end
else
  yScalar = zeros(1,nGroups);
end

n = size(E.Ea{1},1);
chunk = max(1,floor(4e6/n^2));
nOff = numel(E.offsetGrid);
for e0 = 1:chunk:E.nE
  ec = e0:min(e0+chunk-1,E.nE);
  Ed = entrySlice(E,1,ec);
  for iOff = 1:nOff
    if seq.finitePulses
      offsets = Ed.dE0 + E.offsetGrid(iOff)*Ed.lw/(2*sqrt(log(2)));
      PulseSub = finitePulsePropagators(Ed,seq,offsets);
      wOff = E.offsetWeight(iOff);
    else
      PulseSub = [];
      wOff = 1;
    end
    memo = struct('keys',{{}},'vals',{{}});
    for p = 1:nP
      [G,D,BL,BR,memo] = evaluatePathways(Ed,seq,p,PulseSub,memo);
      wp = pref(p)*wOff*reshape(Ed.w,1,1,[]);
      if nDims==0
        % no increments: G is the identity, signal = trace(D)
        v = sum(diag3(D).*reshape(wp,1,[]),'all');
        yScalar(group(p)) = yScalar(group(p)) + v;
      else
        G = G.*wp;
        iG = group(p);
        bRe_ = bRe{iG}; bIm_ = bIm{iG};
        iL = seq.idxIncL(p,:); iR = seq.idxIncR(p,:);
        switch numel(BL)
          case 0
            for k = 1:numel(ec)
              sf_peaks(scheme,bRe_,bIm_,seq.dt,iL,iR,Ed.Ea(:,k),Ed.Eb(:,k),G(:,:,k),D(:,:,k));
            end
          case 1
            L1 = BL{1}; R1 = BR{1};
            for k = 1:numel(ec)
              sf_peaks(scheme,bRe_,bIm_,seq.dt,iL,iR,Ed.Ea(:,k),Ed.Eb(:,k),G(:,:,k),D(:,:,k),...
                L1(:,:,k),R1(:,:,k));
            end
          case 2
            L1 = BL{1}; R1 = BR{1}; L2 = BL{2}; R2 = BR{2};
            for k = 1:numel(ec)
              sf_peaks(scheme,bRe_,bIm_,seq.dt,iL,iR,Ed.Ea(:,k),Ed.Eb(:,k),G(:,:,k),D(:,:,k),...
                L1(:,:,k),L2(:,:,k),R1(:,:,k),R2(:,:,k));
            end
          case 3
            L1 = BL{1}; R1 = BR{1}; L2 = BL{2}; R2 = BR{2}; L3 = BL{3}; R3 = BR{3};
            for k = 1:numel(ec)
              sf_peaks(scheme,bRe_,bIm_,seq.dt,iL,iR,Ed.Ea(:,k),Ed.Eb(:,k),G(:,:,k),D(:,:,k),...
                L1(:,:,k),L2(:,:,k),L3(:,:,k),R1(:,:,k),R2(:,:,k),R3(:,:,k));
            end
        end
      end
    end
  end
end

% Inverse Fourier transform, decays
if nDims==0
  y = sum(yScalar);
  return
end
y = 0;
for iG = 1:nGroups
  buff = complex(bRe{iG},bIm{iG});
  if nDims==1
    td_ = ifft(buff)*numel(buff);
    td_ = td_(1:seq.nPointsInc);
  else
    td_ = ifft(buff,[],1);
    td_ = ifft(td_(1:seq.nPointsInc(1),:),[],2);
    td_ = td_(:,1:seq.nPointsInc(2))*numel(buff);
  end
  y = y + applyDecay(td_,rates(iG,:),seq);
end

end

%-------------------------------------------------------------------------------
function td = applyDecay(td,rate,seq)
if all(rate==0), return; end
if numel(seq.nPointsInc)==1
  td = td.*exp(-rate(1)*seq.dt(1)*(0:seq.nPointsInc-1));
else
  d1 = exp(-rate(1)*seq.dt(1)*(0:seq.nPointsInc(1)-1)).';
  d2 = exp(-rate(2)*seq.dt(2)*(0:seq.nPointsInc(2)-1));
  td = td.*(d1*d2);
end
end

%-------------------------------------------------------------------------------
% Data of entries ec for subspace s, in page form (n x 1 x K, n x n x K)
function Ed = entrySlice(E,s,ec)
K = numel(ec);
n = size(E.Ea{s},1);
Ed.Ea = E.Ea{s}(:,ec);
Ed.Eb = E.Eb{s}(:,ec);
Ed.Ea3 = reshape(Ed.Ea,n,1,K);
Ed.Eb3 = reshape(Ed.Eb,n,1,K);
Ed.M = E.M{s}(:,:,ec);
Ed.Mt = pagectranspose(Ed.M);
Ed.w = E.w(ec);
Ed.dE0 = E.dE0(ec);
Ed.lw = E.lw(ec);
Ed.n = n;
Ed.K = K;
end

%-------------------------------------------------------------------------------
% Pulse sub-propagators of the finite microwave pulses for all entries at
% the given offsets (MHz): blocks 1 = aa, 2 = bb, 3 = ab, 4 = ba
function PulseSub = finitePulsePropagators(Ed,seq,offsets)
N = seq.N;
n = Ed.n; K = Ed.K;
PulseSub = cell(N,4);
idxa = 1:n; idxb = n+1:2*n;
Szdiag = [ones(n,1); -ones(n,1)]/2;
for i = 1:N
  if seq.isRF(i) || seq.Pulse(i).ideal, continue; end
  tp = seq.Pulse(i).tp;
  nu1 = seq.Pulse(i).Flip/(2*pi*tp);  % MHz
  phaseFactor = exp(-1i*seq.Pulse(i).Phase);
  U = zeros(2*n,2*n,K);
  for k = 1:K
    Hp = [diag(Ed.Ea(:,k)), Ed.M(:,:,k)*nu1/2; Ed.Mt(:,:,k)*nu1/2, diag(Ed.Eb(:,k))];
    Hp = Hp + diag((offsets(k)-seq.Pulse(i).Frequency)*Szdiag);
    [V,lam] = eig((Hp+Hp')/2,'vector');
    U(:,:,k) = (V.*exp(-2i*pi*tp*lam.'))*V';
  end
  PulseSub{i,1} = U(idxa,idxa,:);
  PulseSub{i,2} = U(idxb,idxb,:);
  PulseSub{i,3} = U(idxa,idxb,:)*phaseFactor;
  PulseSub{i,4} = U(idxb,idxa,:)*conj(phaseFactor);
end
end

%-------------------------------------------------------------------------------
% Evaluate pathway p for all entries: G (state before the first
% incremented delay), the transfer blocks between incremented delays, and
% D (detection operator S+, propagated back from the end). Operands are
% diagonal (n x 1 x K) or full (n x n x K). Products of prefixes are
% shared between pathways (memo).
function [G,D,BL,BR,memo] = evaluatePathways(Ed,seq,p,PulseSub,memo)
if nargin<5 || isempty(memo), memo = struct('keys',{{}},'vals',{{}}); end
N = seq.N;
incIntervals = find(seq.Inc~=0);
segStart = [1 incIntervals+1];
segEnd = [incIntervals N];
nSeg = numel(segStart);
BL = cell(1,nSeg-2); BR = cell(1,nSeg-2);
G = [];
for iSeg = 1:nSeg
  L = []; R = [];
  keyL = sprintf('L%d:',segStart(iSeg)); keyR = sprintf('R%d:',segStart(iSeg));
  for i = segStart(iSeg):segEnd(iSeg)
    cL = seq.idxPulseL(p,i); cR = seq.idxPulseR(p,i);
    keyL = [keyL char('0'+cL)]; %#ok<AGROW>
    keyR = [keyR char('0'+cR)]; %#ok<AGROW>
    [hit,val] = memoGet(memo,keyL);
    if hit
      L = val;
    else
      L = mulLeft(stepBlock(Ed,seq,PulseSub,i,cL),L);
      memo = memoSet(memo,keyL,L);
    end
    [hit,val] = memoGet(memo,keyR);
    if hit
      R = val;
    else
      R = mulRight(R,stepBlock(Ed,seq,PulseSub,i,cR));
      memo = memoSet(memo,keyR,R);
    end
  end
  if nSeg==1
    % no increments: D = R*M*L
    D = mulFull(mulFull(R,Ed.M,Ed.n),L,Ed.n);
    G = [];
  elseif iSeg==1
    G = mulFull(L,R,Ed.n);
  elseif iSeg==nSeg
    D = mulFull(mulFull(R,Ed.M,Ed.n),L,Ed.n);
  else
    BL{iSeg-1} = toFull(L,Ed.n);
    BR{iSeg-1} = toFull(R,Ed.n);
  end
end
end

function [hit,val] = memoGet(memo,key)
k = find(strcmp(memo.keys,key),1);
hit = ~isempty(k);
if hit, val = memo.vals{k}; else, val = []; end
end
function memo = memoSet(memo,key,val)
memo.keys{end+1} = key;
memo.vals{end+1} = val;
end

% Propagator block of interval i (pulse, then free evolution over the
% initial delay): blocks 1 = aa, 2 = bb, 3 = ab, 4 = ba. The free evolution
% is in the manifold of the row index.
function S = stepBlock(Ed,seq,PulseSub,i,c)
t = seq.t(i);
if seq.isRF(i)
  t = t + seq.Pulse(i).tp;  % free evolution during the RF pulse
end
switch c
  case {1,3}, f = exp(-2i*pi*Ed.Ea3*t);
  case {2,4}, f = exp(-2i*pi*Ed.Eb3*t);
end
if seq.isRF(i) || seq.Pulse(i).ideal
  switch c
    case {1,2}, S = f;
    case 3, S = f.*Ed.M;
    case 4, S = f.*Ed.Mt;
  end
else
  S = f.*PulseSub{i,c};
end
end

% P*L for operands that are empty (identity), diagonal or full
function X = mulLeft(P,L)
if isempty(L), X = P; return; end
if size(P,2)==1 || size(L,2)==1
  if size(L,2)==1 && size(P,2)~=1
    X = P.*permute(L,[2 1 3]);
  else
    X = P.*L;
  end
else
  X = pagemtimes(P,L);
end
end

% R*P' for operands that are empty (identity), diagonal or full
function X = mulRight(R,P)
if size(P,2)==1
  Pt = conj(P);
  if isempty(R), X = Pt; return; end
  if size(R,2)==1
    X = R.*Pt;
  else
    X = R.*permute(Pt,[2 1 3]);
  end
else
  if isempty(R), X = pagectranspose(P); return; end
  if size(R,2)==1
    X = R.*pagectranspose(P);
  else
    X = pagemtimes(R,'none',P,'ctranspose');
  end
end
end

% A*B with diagonal/empty/full operands, full result
function X = mulFull(A,B,n)
if isempty(A) && isempty(B), X = []; return; end
if isempty(A), X = toFull(B,n); return; end
if isempty(B), X = toFull(A,n); return; end
if size(A,2)==1 && size(B,2)==1
  X = toFull(A.*B,n);
elseif size(A,2)==1
  X = A.*B;
elseif size(B,2)==1
  X = A.*permute(B,[2 1 3]);
else
  X = pagemtimes(A,B);
end
end

function X = toFull(A,n)
if isempty(A)
  X = eye(n);
elseif size(A,2)==1 && n>1
  X = eye(n).*permute(A,[2 1 3]);
else
  X = A;
end
end

function d = diag3(X)
n = size(X,1);
K = size(X,3);
X = reshape(X,n*n,K);
d = X(1:n+1:end,:);
end


%===============================================================================
% Symmetry probe for field reversal: evaluates the signal of a random
% small system with M and conj(M) at a few points, without binning, and
% checks whether the reversed signal is the complex conjugate ("conj") or
% the conjugate transpose ("conjtranspose") of the detection-phased signal.
function shortcut = symmetryProbe(seq)

shortcut = "compute";
rs = RandStream('mt19937ar','Seed',1);
n = 3;
[Va,~] = qr(randn(rs,n)+1i*randn(rs,n));
[Vb,~] = qr(randn(rs,n)+1i*randn(rs,n));
Ea = sort(randn(rs,n,1)*3);
Eb = sort(randn(rs,n,1)*3);
nDims = numel(seq.kind1Dims);
if nDims==0, return; end
nTest = 4;
pts = randi(rs,8,nTest,nDims)-1;
S = zeros(nTest,2); Ssw = zeros(nTest,2);
for rev = 1:2
  if rev==1, Ma = Va; Mb = Vb; else, Ma = conj(Va); Mb = conj(Vb); end
  Ed.Ea = Ea; Ed.Eb = Eb; Ed.Ea3 = Ea; Ed.Eb3 = Eb;
  Ed.M = Ma'*Mb; Ed.Mt = Ed.M'; Ed.n = n; Ed.K = 1;
  for it = 1:nTest
    S(it,rev) = pointSignal(Ed,seq,pts(it,:));
    if nDims==2
      Ssw(it,rev) = pointSignal(Ed,seq,pts(it,[2 1]));
    end
  end
end
ph = exp(-1i*seq.DetPhase);
S = S*ph; Ssw = Ssw*ph;
scale = max(abs(S(:)));
if scale==0, return; end
if all(abs(S(:,2)-conj(S(:,1)))<1e-9*scale)
  shortcut = "conj";
elseif nDims==2 && seq.dt(1)==seq.dt(2) && seq.nPointsInc(1)==seq.nPointsInc(2) && ...
    all(abs(Ssw(:,2)-conj(S(:,1)))<1e-9*scale)
  shortcut = "conjtranspose";
end
end

% Signal at one point of the incremented dimensions, without binning
function v = pointSignal(Ed,seq,k)
s = seq;
for i = find(seq.Inc~=0)
  s.t(i) = seq.t(i) + k(seq.Inc(i))*seq.dt(seq.Inc(i));
end
s.Inc(:) = 0;
v = 0;
for p = 1:seq.nPathways
  [~,D] = evaluatePathways(Ed,s,p,[]);
  v = v + seq.pathwayprefactor(p)*trace(D);
end
end


%===============================================================================
% ENDOR: echo as a function of the RF carrier frequency (kind-2 dimension),
% for a single rectangular RF pulse. The RF pulse transfers population
% between the two levels u,v of each NMR transition with probability
% P(Delta) = (nu1/Omega)^2*sin(pi*Omega*tp)^2, Omega = sqrt(nu1^2+Delta^2),
% at the detuning Delta = nu_uv - rfFreq - Frequency. Terms that depend on the absolute RF phase
% average out (RF not phase-locked to the microwave, and the chi average of
% powders), leaving the population transfer.
function y = endorSweep(E,seq,Sys,Exp,NucHams,nSubSpaces)

if nSubSpaces>1
  error('ENDOR is not supported with the product rule.');
end
rfFreqs = seq.rfAxis(:).';
nRfFreqs = numel(rfFreqs);
iRF = find(seq.isRF);
rfPulse = seq.Pulse(iRF);
N = seq.N;
nP = seq.nPathways;

% Sequence split at the RF pulse: intervals 1..iRF-1 (state G before the RF
% pulse) and iRF..N (detection D after it)
if iRF==1
  error('The RF pulse cannot be the first pulse.');
end
% The population-transfer model needs the electron spin in a population
% state during the RF pulse (blocks 1 = aa, 2 = bb, not 3 = ab, 4 = ba)
if any(seq.ctp(:,iRF)>2)
  error(['The subspace engine simulates ENDOR only for pathways with electron coherence order 0 ' ...
    'during the RF pulse. Select them with Exp.Pathways (e.g. ''+00-'' for Mims ENDOR), ' ...
    'or use Opt.Method = ''fullspace''.']);
end
seqG = seq;
seqG.Inc = zeros(1,N);
seqG.Inc(iRF-1) = 1;  % split point: G = state at the end of interval iRF-1

% RF nutation frequency: nominal from the flip angle, or per transition
% from B1 (shaped RF pulses are excluded by sf_parsesequence)
useB1 = ~isempty(rfPulse.B1);
if ~useB1
  nu1Nominal = rfPulse.Flip/(2*pi*rfPulse.tp);  % MHz
end

yBaseline = 0;  % signal without RF pulse
yTransfer = zeros(1,nRfFreqs);  % change of the signal by the RF pulse
n = size(E.Ea{1},1);
chunk = max(1,floor(4e6/n^2));
nOff = numel(E.offsetGrid);
[idxU,idxV] = find(triu(true(n),1));  % levels u<v of each NMR transition
nTransitions = numel(idxU);
for e0 = 1:chunk:E.nE
  ec = e0:min(e0+chunk-1,E.nE);
  Ed = entrySlice(E,1,ec);
  K = numel(ec);
  % RF nutation frequencies (MHz) of the transitions in both manifolds
  if useB1
    nu1 = rfNutation(E,ec,NucHams,Sys,Exp,rfPulse.B1,idxU,idxV); % {a,b} nTransitions x K
  end
  % population transfer amplitudes of the NMR transitions in the manifold
  % of the RF pulse (1 = alpha, 2 = beta), summed over offsets and pathways
  transferAmp = {zeros(nTransitions,K),zeros(nTransitions,K)};
  for iOff = 1:nOff
    if seq.finitePulses
      offsets = Ed.dE0 + E.offsetGrid(iOff)*Ed.lw/(2*sqrt(log(2)));
      PulseSub = finitePulsePropagators(Ed,seq,offsets);
      wOff = E.offsetWeight(iOff);
    else
      PulseSub = []; wOff = 1;
    end
    for p = 1:nP
      [G,D] = evaluatePathways(Ed,seqG,p,PulseSub);
      wp = seq.pathwayprefactor(p)*wOff*reshape(Ed.w,1,1,[]);
      G = G.*wp;
      % baseline: trace(D*G)
      yBaseline = yBaseline + sum(G.*permute(D,[2 1 3]),'all');
      % population transfer
      popG = diag3(G); popD = diag3(D);
      if seq.ctp(p,iRF)==1, iManifold = 1; else, iManifold = 2; end
      transferAmp{iManifold} = transferAmp{iManifold} + ...
        (popG(idxV,:)-popG(idxU,:)).*(popD(idxU,:)-popD(idxV,:));  % nTransitions x K
    end
  end
  % RF response: depends only on the NMR transition frequencies and
  % nutation frequencies, not on the offset or the pathway
  for iManifold = 1:2
    if iManifold==1, Enuc = Ed.Ea; else, Enuc = Ed.Eb; end
    nuNMR = abs(Enuc(idxV,:)-Enuc(idxU,:));
    if useB1
      nu1Trans = nu1{iManifold};
    else
      nu1Trans = nu1Nominal*ones(size(nuNMR));
    end
    amp = transferAmp{iManifold}(:); nuNMR = nuNMR(:); nu1Trans = nu1Trans(:);
    nonzero = abs(amp)>0;
    amp = amp(nonzero); nuNMR = nuNMR(nonzero); nu1Trans = nu1Trans(nonzero);
    % accumulate over the RF axis in chunks
    chunkRF = max(1,floor(2e6/nRfFreqs));
    for j0 = 1:chunkRF:numel(amp)
      idx = j0:min(j0+chunkRF-1,numel(amp));
      Delta = nuNMR(idx) - rfFreqs - rfPulse.Frequency;  % pulse offset from the carrier
      Omega = sqrt(nu1Trans(idx).^2 + Delta.^2);
      Ptransfer = (nu1Trans(idx)./Omega).^2.*sin(pi*Omega*rfPulse.tp).^2;
      yTransfer = yTransfer + amp(idx).'*Ptransfer;
    end
  end
end
y = (seq.rfOn+seq.rfOff)*yBaseline + seq.rfOn*yTransfer;
end

%-------------------------------------------------------------------------------
% RF nutation frequencies 2*B1*|<u|O|v>| (MHz) for the transitions (idxU,idxV)
% of both manifolds, with the effective RF operator O (MHz/mT): nuclear
% Zeeman along the RF field direction plus the hyperfine enhancement from
% the tilt of <S> by the RF field. For powders, |<u|O|v>| is averaged
% (in quadrature) over the RF field direction in the plane perpendicular
% to the static field.
function nu1 = rfNutation(E,ec,NucHams,Sys,Exp,B1,idxU,idxV)
if isempty(E.geff)
  error('B1 RF pulses are only supported for S = 1/2 with all nuclei in the subspace engine.');
end
K = numel(ec);
nTransitions = numel(idxU);
nu1 = {zeros(nTransitions,K),zeros(nTransitions,K)};
nChi = 8;
chiList = (0:nChi-1)*pi/nChi;
g = Sys.g;
if ~Sys.fullg
  Rg = erot(Sys.gFrame).';
  g = Rg*diag(g)*Rg.';
end
for k = 1:K
  e = ec(k);
  geffNorm = E.geffNorm(e);
  q = E.geff(:,e)/geffNorm;  % direction of <S>
  xL = E.xL(:,e); yL = E.yL(:,e);
  for iManifold = 1:2
    if iManifold==1
      V = E.Va{1}(:,:,e); mS = +1/2;
    else
      V = E.Vb{1}(:,:,e); mS = -1/2;
    end
    matElSq = zeros(nTransitions,1);
    for chi = chiList
      nB1 = cos(chi)*xL + sin(chi)*yL;  % RF field direction
      gB1 = g.'*nB1;
      dqdB = (gB1 - q*(q.'*gB1))/(geffNorm*Exp.Field);  % tilt of <S> direction, per mT
      dSdB = mS*dqdB;
      Orf = nB1(1)*NucHams(1).Hnzx + nB1(2)*NucHams(1).Hnzy + nB1(3)*NucHams(1).Hnzz + ...
        dSdB(1)*NucHams(1).Hhfx + dSdB(2)*NucHams(1).Hhfy + dSdB(3)*NucHams(1).Hhfz;
      OrfEig = V'*Orf*V;
      matElSq = matElSq + abs(OrfEig(sub2ind(size(OrfEig),idxU,idxV))).^2;
    end
    nu1{iManifold}(:,k) = 2*B1*sqrt(matElSq/nChi);
  end
end
end
