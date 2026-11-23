% sf_engine_manifold  Manifold engine of saffron
%
%   [y,info] = sf_engine_manifold(Sys,Exp,Opt,Orientations)
%
%   Simulates a pulse EPR signal as a sum over electron coherence transfer
%   pathways, in the nuclear subspaces of the electron spin manifolds. This
%   generalizes the formalism of Mims to arbitrary sequences:
%     W. B. Mims, Phys. Rev. B 5, 2409 (1972), doi:10.1103/PhysRevB.5.2409
%
%   Input:
%     Sys           validated spin system (single isotopologue)
%     Exp           experiment, parsed by saffron; uses
%       Field         magnetic field (mT)
%       mwFreq        microwave carrier frequency (GHz), optional
%       DetFreq       detection frame frequency (GHz), optional, default
%                     mwFreq (needs mwFreq)
%       OriWeights    weights of the orientations, sum 1
%     Opt           options, parsed by saffron; uses
%       fastSeq       parsed sequence, from sf_parsesequence
%       fastExp       canonical experiment, for re-parsing the sequence at
%                     each point of point-by-point dimensions
%       Nuclei, ProductRule, Transitions, OriThreshold, OffsetTolerance,
%       Expand, Binning, Relaxation    see the saffron documentation
%       crystalSample, nOctants   sample type and grid octants
%     Orientations  orientations (phi,theta,chi) in radians, one per row
%
%   Output:
%     y     signal, with one array dimension per indirect dimension
%           (Exp.Dim1, Exp.Dim2, ...), normalized and phased (Exp.DetPhase)
%     info  structure with fields
%       decayRates            decay rates (1/µs) of the pathway groups, for
%                             processing; empty for ENDOR
%       ProductRule           true if the product rule was used
%       reversedField         handling of field-reversed partners: "none",
%                             "compute", "conj" or "conjtranspose"
%       nSkippedOrientations  number of orientations without excitation
%       nEntries              number of entries (see below)
%       binnedDims            dimensions binned in the frequency domain
%                             (the other delay dimensions are stepped)
%       nOffsets              largest number of resonance offsets in the
%                             average over the inhomogeneous line (1: no
%                             line, or the exact average for ideal pulses)
%
%   Physics (see the saffron documentation and the redesign notes):
%   - Each excited EPR transition (b -> a) is a two-level electron system.
%     The nuclear sub-Hamiltonians of the two manifolds use the expectation
%     value <S> of the electron spin (high-field approximation). Their
%     eigenvectors Va and Vb give the overlap matrix M = Va'*Vb.
%   - Initial state: -Sz, or an electron-only state (Sys.initState), with
%     equal populations of the nuclear sublevels. Each transition enters
%     with the population difference of its two levels (in E.w).
%   - Ideal microwave pulses give pathway prefactors (sf_parsesequence);
%     rectangular pulses are exact 2n x 2n propagators at each resonance
%     offset. Free evolution is a diagonal phase factor.
%   - Linearly incremented delays: each interval gets longer by a fixed
%     amount per step (incrementation scheme, seq.incScheme). Along a
%     binned dimension, the signal
%     of each pathway is a sum of exponentials, binned in the frequency
%     domain (sf_peaks) and Fourier transformed. Along a stepped dimension,
%     the delays are fixed at each point and the signal is computed point
%     by point (Opt.Binning, chooseBinning).
%   - Inhomogeneous line (Sys.lwAnisoG): a Gaussian distribution of
%     resonance offsets per entry, with 1/e half width sigma. With ideal
%     pulses, the offset enters a pathway only as a phase over its net
%     electron coherence time Tnet, and the average over the line is
%     exp(-(pi*sigma*Tnet)^2). This is used wherever the signal is not
%     binned (single points, product rule, ENDOR). Otherwise, the average
%     uses the trapezoid rule (offsetQuadrature).
%   - RF pulses (ENDOR, rectangular only): population transfer between the
%     levels of each NMR transition, with the two-level response (Rabi
%     formula) at the detuning from the RF carrier, averaged over the RF
%     phase. The nutation frequency of each transition follows from its
%     matrix element of the RF operator (nuclear Zeeman interaction and,
%     by default, the hyperfine enhancement, Opt.rfEnhancement) and the
%     RF field B2 (Pulse.B2, or nominal from Pulse.Flip). In powders, the
%     response is averaged over the RF field direction in the plane
%     perpendicular to the static field (see endorSweep).
%   - Relaxation: exp(-t/T2) for intervals with electron coherence and
%     exp(-t/T1) for intervals with electron populations, over the delays
%     with relaxation on (Opt.Relaxation, one value or one per element of
%     Exp.Sequence); none during pulses.
%
%   Data structures of the local functions:
%     E  entries, one per excited EPR transition and orientation (and
%        reversed field), with n nuclear states and nE entries
%       nE                number of entries
%       w                 weight of each entry (orientation weight, selection
%                         weight, population difference of the initial
%                         state)
%       dE0               EPR resonance offset from Exp.mwFreq (MHz)
%       lw                inhomogeneous FWHM (MHz)
%       ori               orientation index
%       rev               true for field-reversed partners (the second half
%                         of the entries, if present)
%       moment            transition moments |<a|S+|b>|, empty: all 1
%       Ea{s}, Eb{s}      nuclear energies (MHz) of manifolds a and b, one
%                         cell per subspace s, n x nE
%       M{s}              overlap matrices Va'*Vb, n x n x nE
%       Va{s}, Vb{s}      nuclear eigenvectors, n x n x nE (ENDOR only)
%       xL, yL, zL        lab axes in the molecular frame, 3 x nE (ENDOR
%                         only; zL only for the general electron system)
%       dSa, dSb          tilt of the electron spin <S> (molecular frame)
%                         in manifolds a and b by an RF field (per mT)
%                         along xL (rows 1-3) and yL (rows 4-6), 6 x nE;
%                         zero without hyperfine enhancement (ENDOR only)
%       offCenter, offWidth   center and FWHM (MHz) of the resonance offset
%                         distribution of the excited spin packets
%       offsetTolerance   error target of the average over this
%                         distribution (Opt.OffsetTolerance), see
%                         offsetQuadrature
%     Ed  data of K entries of one subspace, from entrySlice
%       Ea, Eb            nuclear energies (MHz), n x K
%       Ea3, Eb3          same, n x 1 x K
%       M, Mt             overlap matrices and their conjugate transposes,
%                         n x n x K
%       w, dE0, offCenter, offWidth, moment   as in E, 1 x K
%       n, K              number of nuclear states and of entries
%       original          position of the original entry of a field-reversed
%                         partner in this slice, 0 if not present, 1 x K

function [y,info] = sf_engine_manifold(Sys,Exp,Opt,Orientations)

seq = Opt.fastSeq;
nOrientations = size(Orientations,1);
info = struct;

% Relaxation times (0 or not given: no relaxation)
T1 = 0; T2 = 0;
if isfield(Sys,'T1') && ~isempty(Sys.T1)
  if numel(Sys.T1)~=1
    error('The manifold engine does not support transition-selective relaxation times (Sys.T1). Use Opt.Method = ''full''.');
  end
  T1 = Sys.T1;
end
if isfield(Sys,'T2') && ~isempty(Sys.T2)
  if numel(Sys.T2)~=1
    error('The manifold engine does not support transition-selective relaxation times (Sys.T2). Use Opt.Method = ''full''.');
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
rfEnhancement = ~isfield(Opt,'rfEnhancement') || Opt.rfEnhancement;

% Nuclei included in the nuclear subspaces
%-------------------------------------------------------------------------------
% Nuclear-nuclear couplings (Sys.nn): coupled pairs and their tensors in
% the molecular frame. Coupled nuclei must be in the nuclear subspace,
% since the core system cannot hold these couplings.
nnPairs = zeros(0,2);
nnTensors = {};
if Sys.nNuclei>=2 && any(Sys.nn(:))
  allPairs = nchoosek(1:Sys.nNuclei,2);
  for iPair = 1:size(allPairs,1)
    if Sys.fullnn
      J = Sys.nn(3*(iPair-1)+(1:3),:);
    else
      J = diag(Sys.nn(iPair,:));
    end
    if ~any(J(:)), continue; end
    R_nn2M = erot(Sys.nnFrame(iPair,:)).';  % nn frame -> molecular frame
    nnPairs(end+1,:) = allPairs(iPair,:); %#ok<AGROW>
    nnTensors{end+1} = R_nn2M*J*R_nn2M.'; %#ok<AGROW>
  end
end
nnCoupled = ismember(1:Sys.nNuclei,nnPairs(:));
if isfield(Opt,'Nuclei') && ~isempty(Opt.Nuclei)
  if any(Opt.Nuclei>Sys.nNuclei) || any(Opt.Nuclei<1)
    error('Index in Opt.Nuclei out of range.');
  end
  shfNuclei = Opt.Nuclei;
  missing = find(nnCoupled & ~ismember(1:Sys.nNuclei,shfNuclei));
  if ~isempty(missing)
    error('Opt.Nuclei must include all nuclei with nuclear-nuclear couplings (Sys.nn). Missing: %s.',...
      num2str(missing));
  end
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
    if any(strong & nnCoupled)
      error(['Nucleus %s has a hyperfine coupling larger than the excitation window of Opt.IdealPulses ' ...
        'and a nuclear-nuclear coupling (Sys.nn), so it cannot be moved to the core system. ' ...
        'Use finite pulses (Opt.IdealPulses = false).'],num2str(find(strong & nnCoupled)));
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
  if productRule && any(nnCoupled)
    error('Opt.ProductRule requires independent nuclei and cannot be used with nuclear-nuclear couplings (Sys.nn).');
  end
end
if numel(shfNuclei)<=1, productRule = false; end
info.ProductRule = productRule;

% Nuclear sub-Hamiltonians
%-------------------------------------------------------------------------------
NucHams = nuclearHamiltonians(Sys,shfNuclei,productRule,nnPairs,nnTensors);
nSubSpaces = numel(NucHams);

% Electronic core system (not needed for S = 1/2 with all nuclei)
if ~twoElectronManifolds
  % all nn-coupled nuclei are in the nuclear subspace, so the core system
  % has no nuclear-nuclear couplings
  coreSys = rmfield(Sys,intersect(fieldnames(Sys),{'nn','nnFrame'}));
  coreSys = nucspinrmv(coreSys,shfNuclei);
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

% Initial state (Sys.initState): an electron-only state, expanded with equal
% populations of the nuclear sublevels (sf_selectengine rejects states with
% a nuclear part). Only the populations of the two levels of each EPR
% transition enter, as the weight rho_bb - rho_aa of the transition; for
% -Sz, this is <a|Sz|a> - <b|Sz|b>. These populations of the eigenstates
% do not change when the molecule is rotated about the field (third Euler
% angle chi), while the coherences rho_ab between the two levels average
% out over chi in powders (sf_selectengine rejects crystals). The
% full engine instead removes all coherences between different
% electron mS from the state, which is the same with its secular
% Hamiltonian, whose eigenstates have pure mS along the field. The two
% differ by the deviation of the eigenstates from pure mS: to first order
% in D/nu_e for S > 1/2, and for S = 1/2 with anisotropic g by the factor
% q.zL, the cosine of the angle between the quantization axis q and the
% field. Eigenstates ('eigen' basis) are those of the electron
% Hamiltonian without hyperfine couplings, to first order as in the
% nuclear sub-Hamiltonians. The states are in the molecular frame, with
% the trace of the electron state divided by the number of nuclear
% states, as in the full engine.
userInitState = ~isempty(Sys.initState);
if userInitState
  SysE = rmfield(Sys,intersect(fieldnames(Sys),{'nn','nnFrame'}));
  SysE = nucspinrmv(SysE,1:Sys.nNuclei);
  SysE.processed = false;
  [SysE,err] = validatespinsys(SysE);
  error(err);
  nNucStates = prod(2*Sys.I+1);
  % 'eigen': the state depends on the field direction, through the
  % eigenstates of the electron Hamiltonian H0e - B*mu.zL; otherwise it is
  % the same state rhoE at all orientations
  eigenBasis = strcmp(Sys.initState{2},'eigen');
  if eigenBasis
    [H0e,muxE,muyE,muzE] = ham(SysE);
  else
    rhoE = initstate(SysE)/nNucStates;
  end
end

OrientationSelection = isfield(Exp,'mwFreq') && ~isempty(seq.ExciteWidth);
if ~isfield(Opt,'OriThreshold'), Opt.OriThreshold = 0.005; end
if ~isfield(Sys,'lwAnisoG'), Sys.lwAnisoG = [0 0 0]; end
lwAnisoG = Sys.lwAnisoG(:).'.*[1 1 1];

% Field-reversed orientations, for partial-sphere powder grids. Reversing
% the field direction keeps the nuclear energies and conjugates the overlap
% matrices, M -> conj(M). Without static nuclear terms (Q, nn), M can be
% made real for each nucleus, and the reversed orientations give the same
% signal. With them, the reversed partners are needed (see
% fieldReversalShortcut below). The shortcut is checked on the parsed
% sequence, i.e. with the pulses and delays of Exp.Sequence as given,
% without the changes of the point-by-point (kind 3) dimensions. These can
% break the symmetry (e.g. a changed delay), so with them the partners are
% computed.
% An initial state can also differ between n and -n (time-odd states, e.g.
% an electron polarization along a molecular axis): the weights of the
% reversed orientations are computed (reversedWeights). The shortcuts need
% equal weights, so with an initial state they are chosen once the weights
% are known (below). The full engine treats the same cases (fieldReversal
% in sf_engine_full).
partialSphereGrid = ~Opt.crystalSample && isfield(Opt,'nOctants') && Opt.nOctants~=8;
hasStaticTerms = nSubSpaces>0 && any(arrayfun(@(N)any(N.Hstatic(:)~=0),NucHams));
reversedWeights = partialSphereGrid && userInitState;
shortcutAllowed = ~seq.finitePulses && ~isENDOR && ~any(seq.kind==3);
reversedField = "none";
if partialSphereGrid && hasStaticTerms
  reversedField = "compute";
  if shortcutAllowed && ~userInitState
    reversedField = fieldReversalShortcut(seq);
  end
end

% Entries: one per excited EPR transition and orientation (and reversed
% field), with nuclear energies and overlap matrices of all subspaces
%-------------------------------------------------------------------------------
E = struct;
E.w = zeros(1,0);            % weight (orientation, selection, population difference)
E.dE0 = zeros(1,0);          % EPR resonance offset from mwFreq (MHz)
E.lw = zeros(1,0);           % inhomogeneous FWHM (MHz)
E.ori = zeros(1,0);
E.rev = false(1,0);
E.moment = [];               % transition moments |<a|S+|b>|, empty: all 1 (S = 1/2)
for s = 1:nSubSpaces
  E.Ea{s} = []; E.Eb{s} = []; E.M{s} = [];
  if isENDOR, E.Va{s} = []; E.Vb{s} = []; end
end
if twoElectronManifolds
  zL = zeros(3,nOrientations); xL = zL; yL = zL;
  for iOri = 1:nOrientations
    [xL(:,iOri),yL(:,iOri),zL(:,iOri)] = erot(Orientations(iOri,:),'rows');
  end
  % Zeeman term bmagn*B*zL.'*g*S: the electron spin is quantized along
  % g.'*zL (also for a non-symmetric g)
  geff = g.'*zL;
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
  if userInitState
    % levels a and b: electron spin along +q and -q (q = geff/|geff|), so
    % that rho_bb - rho_aa = -2*Re(trace(rho*(q.S))); the reversed
    % orientation -n has -q. The eigenstates of the 'eigen' basis are the
    % two levels themselves, so their weight is the same at all
    % orientations, and at -n.
    Sops = {full(sop(1/2,'x')),full(sop(1/2,'y')),full(sop(1/2,'z'))};
    wRev = E.w;
    for k = 1:nE
      iOri = oriList(k);
      q = geff(:,iOri)/geffNorm(iOri);
      qS = q(1)*Sops{1} + q(2)*Sops{2} + q(3)*Sops{3};
      if eigenBasis
        if k==1
          zL1 = zL(:,iOri);
          rho = initstate(SysE,H0e-Exp.Field*(zL1(1)*muxE+zL1(2)*muyE+zL1(3)*muzE))/nNucStates;
          dPop = -2*real(trace(rho*qS));
        end
        E.w(k) = E.w(k)*dPop;
        wRev(k) = wRev(k)*dPop;
      else
        E.w(k) = E.w(k)*(-2)*real(trace(rhoE*qS));
        wRev(k) = wRev(k)*2*real(trace(rhoE*qS));
      end
    end
  end
  if isfield(Exp,'mwFreq')
    E.dE0 = nuEPR(oriList) - 1000*Exp.mwFreq;
  else
    E.dE0 = zeros(1,nE);
  end
  E.lw = lwInh(oriList);
  E.ori = oriList;
  E.rev = false(1,nE);
  for s = 1:nSubSpaces
    n = size(NucHams(s).Hnzz,1);
    Ea = zeros(n,nE); Eb = zeros(n,nE); M = zeros(n,n,nE);
    if isENDOR, VaAll = M; VbAll = M; end
    for k = 1:nE
      iOri = oriList(k);
      q = geff(:,iOri).'/geffNorm(iOri);
      Hnuc = Exp.Field*(zL(1,iOri)*NucHams(s).Hnzx + zL(2,iOri)*NucHams(s).Hnzy + zL(3,iOri)*NucHams(s).Hnzz) + NucHams(s).Hstatic;
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
    E.xL = xL(:,oriList); E.yL = yL(:,oriList);
    % Tilt of <S> (+q/2 in manifold a, -q/2 in b) by an RF field (per mT)
    % along xL (rows 1-3) and yL (rows 4-6), for the hyperfine enhancement
    % (see rfNutation): the quantization axis q = geff/|geff| turns toward
    % g.'*n (n: RF field direction), first order in 1/nu_e (see sf_spintilt)
    E.dSa = zeros(6,nE);
    if rfEnhancement
      q = geff(:,oriList)./geffNorm(oriList);
      gx = g.'*E.xL; gy = g.'*E.yL;
      dqx = (gx - q.*sum(q.*gx,1))./(geffNorm(oriList)*Exp.Field);
      dqy = (gy - q.*sum(q.*gy,1))./(geffNorm(oriList)*Exp.Field);
      E.dSa = [dqx; dqy]/2;
    end
    E.dSb = -E.dSa;
  end
else
  % General electron system: transitions from the core Hamiltonian. Each
  % transition b -> a (a above b) is a two-level system. Its microwave
  % coupling is the co-rotating field component, as in the full
  % engine (rotating-wave approximation): the nutation frequency is
  % nu1*|<a|S+|b>|, with S+ in the lab frame. The initial state -Sz gives
  % the population difference <a|Sz|a> - <b|Sz|b>. Ideal pulses (only
  % S = 1/2 here, with nuclei in the core system) use the nominal flip
  % angle, since |<a|S+|b>| = 1 for the allowed transitions up to the
  % hyperfine mixing in the core system.
  % With an initial state, the core system holds the state of the electron
  % and the core nuclei. The eigenstates at the reversed orientation -n are
  % the time-reversed ones, Theta*|a>, with Theta = Uy*K (K: complex
  % conjugation, Uy = exp(-1i*pi*Sy) of all spins), since the zero-field
  % Hamiltonian is time-even and the Zeeman terms are time-odd. Their
  % populations of the state rho at -n are those of conj(Uy'*rho*Uy) in
  % the eigenstates at n.
  if userInitState
    nCoreNuc = size(H0,1)/prod(2*Sys.S+1);
    if reversedWeights
      SyTotal = 0;
      for iSpin = 1:1+Sys.nNuclei-numel(shfNuclei)
        SyTotal = SyTotal + sop(coreSys,[iSpin,2]);
      end
      Uy = expm(-1i*pi*full(SyTotal));
    end
    if ~eigenBasis
      % the same state at all orientations
      rhoCore = kron(rhoE,eye(nCoreNuc));
      if reversedWeights, rhoCoreRev = conj(Uy'*rhoCore*Uy); end
    end
  end
  rows = {};
  wRev = zeros(1,0);  % weights at -n (reversedWeights)
  for iOri = 1:nOrientations
    [xLab_M,yLab_M,zLab_M] = erot(Orientations(iOri,:),'rows');
    muzL = zLab_M(1)*mux + zLab_M(2)*muy + zLab_M(3)*muz;
    H = H0 - Exp.Field*muzL;
    [eV,eE] = eig(H);
    eE = real(diag(eE));
    SpLab = (xLab_M(1)+1i*yLab_M(1))*Sx + (xLab_M(2)+1i*yLab_M(2))*Sy + (xLab_M(3)+1i*yLab_M(3))*Sz;
    transMoment = abs(eV'*SpLab*eV);  % |<a|S+|b>|
    selWeight = transMoment;
    maxMoment = max(transMoment(:));
    lwInh2 = sum((lwAnisoG(:).*zLab_M(:)).^2);
    if OrientationSelection
      lw2 = lwInh2 + seq.ExciteWidth^2;
      dE = eE - eE.' - Exp.mwFreq*1e3;
      selWeight = selWeight.*exp(-4*log(2)*dE.^2/lw2)/sqrt(1+lwInh2/seq.ExciteWidth^2);
    end
    selWeight(selWeight<Opt.OriThreshold*maxMoment) = 0;
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
      selWeight(rmvTransition) = 0;
    end
    [v,u,wSel] = find(tril(selWeight,-1));
    % Tilt of the electron spin by an RF field (per mT) along xL and yL, for
    % the hyperfine enhancement of RF pulses (see rfNutation)
    if isENDOR && rfEnhancement && ~isempty(u)
      Szf = zLab_M(1)*Sx + zLab_M(2)*Sy + zLab_M(3)*Sz;
      dSx = sf_spintilt(H,-(xLab_M(1)*mux + xLab_M(2)*muy + xLab_M(3)*muz),{Sx,Sy,Sz},Szf);
      dSy = sf_spintilt(H,-(yLab_M(1)*mux + yLab_M(2)*muy + yLab_M(3)*muz),{Sx,Sy,Sz},Szf);
      dSops = [dSx dSy];
    end
    if userInitState && eigenBasis && ~isempty(u)
      % initial state in the core system, from the eigenstates of the
      % electron Hamiltonian at this orientation and, time-reversed, at the
      % reversed one (field -zL)
      muzE_L = zLab_M(1)*muxE + zLab_M(2)*muyE + zLab_M(3)*muzE;
      rhoCore = kron(initstate(SysE,H0e-Exp.Field*muzE_L)/nNucStates,eye(nCoreNuc));
      if reversedWeights
        rhoCoreRev = kron(initstate(SysE,H0e+Exp.Field*muzE_L)/nNucStates,eye(nCoreNuc));
        rhoCoreRev = conj(Uy'*rhoCoreRev*Uy);
      end
    end
    for iT = 1:numel(u)
      b = u(iT); a = v(iT);
      vecA = eV(:,a); vecB = eV(:,b);
      SA = real([vecA'*Sx*vecA; vecA'*Sy*vecA; vecA'*Sz*vecA]);
      SB = real([vecB'*Sx*vecB; vecB'*Sy*vecB; vecB'*Sz*vecB]);
      if userInitState
        dPop = real(vecB'*rhoCore*vecB - vecA'*rhoCore*vecA);
        if reversedWeights
          dPopRev = real(vecB'*rhoCoreRev*vecB - vecA'*rhoCoreRev*vecA);
          wRev(end+1) = Exp.OriWeights(iOri)*wSel(iT)*dPopRev; %#ok<AGROW>
        end
      else
        dPop = zLab_M(:).'*(SA-SB);  % population difference from -Sz
      end
      dSA = zeros(6,1); dSB = zeros(6,1);
      if isENDOR && rfEnhancement
        for c = 1:6
          dSA(c) = real(vecA'*dSops{c}*vecA);
          dSB(c) = real(vecB'*dSops{c}*vecB);
        end
      end
      rows(end+1,:) = {iOri, Exp.OriWeights(iOri)*wSel(iT)*dPop, eE(a)-eE(b), ...
        sqrt(lwInh2), SA, SB, zLab_M, transMoment(a,b), dSA, dSB, xLab_M, yLab_M}; %#ok<AGROW>
    end
  end
  nE = size(rows,1);
  info.nSkippedOrientations = nOrientations - numel(unique([rows{:,1}]));
  E.w = [rows{:,2}];
  E.moment = [rows{:,8}];  % transition moments |<a|S+|b>|
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
      Hnuc = Exp.Field*(zLab_M(1)*NucHams(s).Hnzx + zLab_M(2)*NucHams(s).Hnzy + zLab_M(3)*NucHams(s).Hnzz) + NucHams(s).Hstatic;
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
    E.xL = [rows{:,11}]; E.yL = [rows{:,12}];
    E.dSa = [rows{:,9}]; E.dSb = [rows{:,10}];
  end
end

% Field-reversed partners, computed explicitly (M, Va, Vb -> conj). With
% an initial state, the reversed orientations have their own weights
% wRev. Without partners (no static nuclear terms), both orientations give
% the same signal per unit weight, so the weights are averaged. With equal
% weights (time-even states, e.g. populations of eigenstates), the
% shortcuts apply as without initial state.
if ~reversedWeights, wRev = E.w; end
wScale = max([abs(E.w) abs(wRev) 0]);
if reversedWeights && reversedField=="none"
  E.w = (E.w + wRev)/2;
  wRev = E.w;
elseif reversedWeights && shortcutAllowed && all(abs(wRev-E.w)<=1e-12*wScale)
  reversedField = fieldReversalShortcut(seq);
end
info.reversedField = reversedField;

% Entries without weight do not contribute: transitions between levels with
% equal populations of the initial state, and orientations whose weights at
% n and -n cancel (time-odd states, without partners). At least one entry
% is kept.
if userInitState && nE>0
  keep = abs(E.w)>1e-12*wScale | abs(wRev)>1e-12*wScale;
  if ~any(keep), keep(1) = true; end
  if ~all(keep)
    rowFields = intersect(fieldnames(E),{'w','dE0','lw','ori','rev','moment','dSa','dSb','xL','yL','zL'});
    for iField = 1:numel(rowFields)
      fld = rowFields{iField};
      if ~isempty(E.(fld)), E.(fld) = E.(fld)(:,keep); end
    end
    wRev = wRev(keep);
    for s = 1:nSubSpaces
      E.Ea{s} = E.Ea{s}(:,keep);
      E.Eb{s} = E.Eb{s}(:,keep);
      E.M{s} = E.M{s}(:,:,keep);
      if isENDOR
        E.Va{s} = E.Va{s}(:,:,keep);
        E.Vb{s} = E.Vb{s}(:,:,keep);
      end
    end
    nE = nnz(keep);
  end
end

if reversedField=="compute"
  E.w = [E.w wRev]/2;
  E.dE0 = [E.dE0 E.dE0];
  E.lw = [E.lw E.lw];
  E.ori = [E.ori E.ori];
  E.rev = [false(1,nE) true(1,nE)];
  if ~isempty(E.moment), E.moment = [E.moment E.moment]; end
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
    oriFields = intersect(fieldnames(E),{'dSa','dSb','xL','yL','zL'});
    for iField = 1:numel(oriFields)
      fld = oriFields{iField};
      if ~isempty(E.(fld)), E.(fld) = [E.(fld) E.(fld)]; end
    end
  end
  nE = 2*nE;
end
E.nE = nE;
info.nEntries = nE;

% Resonance offsets of the excited spin packets of each entry: a Gaussian
% distribution with FWHM from Sys.lwAnisoG, around the resonance offset
% dE0. With ideal pulses and an excitation window of FWHM W
% (Opt.IdealPulses), the excited packets follow the product of line and
% window, with center dE0*W^2/(lw^2+W^2) and FWHM lw*W/sqrt(lw^2+W^2) (the
% window weight is in E.w). With finite pulses, the pulses select the
% packets from the whole line.
if ~seq.finitePulses && OrientationSelection
  W2 = seq.ExciteWidth^2;
  E.offCenter = E.dE0.*W2./(E.lw.^2+W2);
  E.offWidth = E.lw*seq.ExciteWidth./sqrt(E.lw.^2+W2);
else
  E.offCenter = E.dE0;
  E.offWidth = E.lw;
end
% The average over the distribution is computed for each parsed sequence,
% with the number of offsets that its timing needs (offsetQuadrature)
E.offsetTolerance = Opt.OffsetTolerance;

% Simulation
%-------------------------------------------------------------------------------
Rel.T1 = T1; Rel.T2 = T2; Rel.on = relaxation;

% Linear delay increments (kind 1): binned in the frequency domain or
% stepped point by point, per dimension (Opt.Binning)
if ~isfield(Opt,'Binning') || isempty(Opt.Binning), Opt.Binning = 'auto'; end
nStates = max([1 cellfun(@(Ea)size(Ea,1),E.Ea)]);
binnedDims = chooseBinning(Opt.Binning,seq,nStates,E.nE,productRule);
seq.binned = binnedDims(seq.kind1Dims);
info.binnedDims = find(binnedDims);
expandFactor = 6;  % 2^expandFactor buffer expansion for binning
if sum(seq.binned)==2, expandFactor = 2; end
if isfield(Opt,'Expand') && ~isempty(Opt.Expand), expandFactor = Opt.Expand; end

if isENDOR
  % RF field amplitude B2 (mT): given, or nominal from the flip angle for
  % the reference nucleus (rfGammaRef, see saffron)
  rfPulse = seq.Pulse(seq.isRF);
  if ~isempty(rfPulse.B2)
    rf.B2 = rfPulse.B2;
  elseif isfield(Opt,'rfGammaRef')
    rf.B2 = rfPulse.Flip/(2*pi*rfPulse.tp*Opt.rfGammaRef);
  else
    error('RF pulses need B2 (in mT), or a nucleus to define a nominal flip angle.');
  end
  rf.tp = rfPulse.tp;
  rf.averageDirections = ~Opt.crystalSample;
  [y,nOffsets,info.nRFDirections] = endorSweep(E,seq,NucHams,nSubSpaces,rf);
  rates = [];
else
  kind3 = find(seq.kind==3);
  if isempty(kind3)
    [y,rates,nOffsets] = pathwaySignal(E,seq,Rel,productRule,expandFactor);
  else
    % Point-by-point dimensions: re-parse the sequence with modified pulses
    % and delays for each point
    nPts = seq.nPoints(kind3);
    yk = cell([nPts 1]);
    nOffsets = 1;
    for iPt = 1:prod(nPts)
      sub = cell(1,numel(kind3));
      [sub{:}] = ind2sub([nPts 1],iPt);
      Expk = Opt.fastExp;
      if ~isempty(seq.echoOrders), Expk.Pathways = seq.echoOrders; end
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
      idealPulses = isfield(Opt,'IdealPulses') && Opt.IdealPulses;
      [seqk,reasons] = sf_parsesequence(Expk,idealPulses,Opt.Relaxation);
      if ~isempty(reasons)
        error('The manifold engine cannot simulate point %d: %s.',iPt,strjoin(reasons,', '));
      end
      seqk.binned = binnedDims(keepDims(seqk.kind1Dims));
      [yk{iPt},rates,nOffk] = pathwaySignal(E,seqk,Rel,productRule,expandFactor);
      nOffsets = max(nOffsets,nOffk);
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
info.nOffsets = nOffsets;
if nOffsets>1
  logmsg(1,'  average over the inhomogeneous line: %d resonance offsets (Opt.OffsetTolerance = %g)',...
    nOffsets,E.offsetTolerance);
end

% Detection phase, field-reversed partners, normalization
y = y*exp(-1i*seq.DetPhase);
switch reversedField
  case "conj"
    y = real(y);
  case "conjtranspose"
    y = (y + y')/2;
end

% Detection frame (Exp.DetFreq): the signal is computed in the frame rotating
% at Exp.mwFreq. In a frame rotating at DetFreq, the detected -1 coherence
% has the phase factor exp(-2i*pi*(DetFreq-mwFreq)*t), with t the time from
% the start of the sequence to the detection point. The factor is the same
% for all orientations and field-reversed partners.
if isfield(Exp,'DetFreq') && ~isempty(Exp.DetFreq) && Exp.DetFreq~=Exp.mwFreq
  t = seq.tLength;
  for d = 1:numel(seq.nPoints)
    t = t + reshape(seq.tLengthChange{d},[ones(1,d-1) seq.nPoints(d) 1]);
  end
  y = y.*reshape(exp(-2i*pi*(Exp.DetFreq-Exp.mwFreq)*1e3*t),size(y));  % GHz*µs -> rad
end
% Normalization: echo amplitude Tr(S+*rho)/Tr(S+'*S+) for the initial state
% -Sz, as in the full engine. Each transition contributes with its
% population difference times its transition moment (in E.w), and
% Tr(S+'*S+) = prod(2I+1)*(2/3)*S*(S+1)*(2S+1), which is prod(2I+1) for S = 1/2.
S = Sys.S;
y = y/(prod(2*Sys.I+1)*(4*S*(S+1)*(2*S+1))/6);

end


%===============================================================================
% nuclearHamiltonians  Build the nuclear sub-Hamiltonians
%
%   NucHams = nuclearHamiltonians(Sys,shfNuclei,productRule,nnPairs,nnTensors)
%
%   Without nuclei in the subspace (no nuclei in the spin system, or all of
%   them in the core system), the nuclear subspace has a single state: each
%   EPR transition is a plain two-level system.
%
%   Input:
%     Sys          validated spin system
%     shfNuclei    indices of the nuclei in the nuclear subspaces
%     productRule  true: one subspace per nucleus; false: one subspace with
%                  all nuclei
%     nnPairs      nucleus indices of the coupled pairs, nPairs x 2
%     nnTensors    molecular-frame nn tensor (3x3, MHz) of each pair, cell
%                  array
%
%   Output:
%     NucHams  structure array with one element per subspace, with fields
%       Hnzx, Hnzy, Hnzz  nuclear Zeeman terms (MHz/mT) along the molecular
%                         axes, with chemical shielding (sigma) and gnscale,
%                         as in ham_nz
%       Hhfx, Hhfy, Hhfz  hyperfine terms (MHz) per unit of the electron
%                         spin expectation value along the molecular axes
%       Hstatic           terms (MHz) that depend neither on the field nor
%                         on the electron state: nuclear quadrupole (Q) and
%                         nuclear-nuclear couplings (nn)
function NucHams = nuclearHamiltonians(Sys,shfNuclei,productRule,nnPairs,nnTensors)

NucHams = struct('Hnzx',{},'Hnzy',{},'Hnzz',{},'Hhfx',{},'Hhfy',{},'Hhfz',{},'Hstatic',{});
I = Sys.I(shfNuclei);
if ~productRule
  dim = prod(2*I+1);
  z = zeros(dim);
  NucHams(1) = struct('Hnzx',z,'Hnzy',z,'Hnzz',z,'Hhfx',z,'Hhfy',z,'Hhfz',z,'Hstatic',z);
end
for ii = 1:numel(shfNuclei)
  iNuc = shfNuclei(ii);
  if productRule
    Ix = sop(I(ii),[1,1]); Iy = sop(I(ii),[1,2]); Iz = sop(I(ii),[1,3]);
    s = ii;
    z = zeros(size(Ix));
    NucHams(s) = struct('Hnzx',z,'Hnzy',z,'Hnzz',z,'Hhfx',z,'Hhfy',z,'Hhfz',z,'Hstatic',z);
  else
    Ix = sop(I,[ii,1]); Iy = sop(I,[ii,2]); Iz = sop(I,[ii,3]);
    s = 1;
  end
  Ix = full(Ix); Iy = full(Iy); Iz = full(Iz);
  % Nuclear Zeeman interaction -gn*muN*B*sigma*I, with the chemical
  % shielding tensor sigma (1 = no shielding) in the molecular frame
  if Sys.fullsigma
    sigma = Sys.sigma((iNuc-1)*3+(1:3),:);
  else
    R_CS2M = erot(Sys.sigmaFrame(iNuc,:)).';
    sigma = R_CS2M*diag(Sys.sigma(iNuc,:))*R_CS2M.';
  end
  pre = -Sys.gn(iNuc)*Sys.gnscale(iNuc)*nmagn/1e3/planck/1e6;  % MHz/mT
  NucHams(s).Hnzx = NucHams(s).Hnzx + pre*(sigma(1,1)*Ix + sigma(1,2)*Iy + sigma(1,3)*Iz);
  NucHams(s).Hnzy = NucHams(s).Hnzy + pre*(sigma(2,1)*Ix + sigma(2,2)*Iy + sigma(2,3)*Iz);
  NucHams(s).Hnzz = NucHams(s).Hnzz + pre*(sigma(3,1)*Ix + sigma(3,2)*Iy + sigma(3,3)*Iz);
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
    NucHams(s).Hstatic = NucHams(s).Hstatic + ...
      Ix*(Q(1,1)*Ix + Q(1,2)*Iy + Q(1,3)*Iz) + ...
      Iy*(Q(2,1)*Ix + Q(2,2)*Iy + Q(2,3)*Iz) + ...
      Iz*(Q(3,1)*Ix + Q(3,2)*Iy + Q(3,3)*Iz);
  end
end

% Nuclear-nuclear couplings I1*nn*I2, the same in both electron manifolds.
% The caller ensures that all coupled nuclei are in shfNuclei and that the
% product rule is off.
for iPair = 1:size(nnPairs,1)
  J = nnTensors{iPair};
  i1 = find(shfNuclei==nnPairs(iPair,1));
  i2 = find(shfNuclei==nnPairs(iPair,2));
  for c1 = 1:3
    I1 = full(sop(I,[i1,c1]));
    for c2 = 1:3
      NucHams(1).Hstatic = NucHams(1).Hstatic + J(c1,c2)*I1*full(sop(I,[i2,c2]));
    end
  end
end
end


%===============================================================================
% pathwaySignal  Signal over the linearly incremented dimensions (kind 1)
%
%   [y,rates,nOff] = pathwaySignal(E,seq,Rel,productRule,expandFactor)
%
%   Sums the signal over all entries, resonance offsets and pathways. Along
%   the binned dimensions (seq.binned), the peaks of each pathway are binned
%   in the frequency domain (sf_peaks) and Fourier transformed. Along the
%   stepped dimensions, the signal is computed point by point: at each
%   point, the incremented delays have fixed lengths (steppedSignal). With
%   ideal pulses, the average over the inhomogeneous line is exact where
%   the signal is not binned (stepped points, product rule), see the header
%   of sf_engine_manifold.
%
%   Input:
%     E             entries, see the header of sf_engine_manifold
%     seq           parsed sequence, from sf_parsesequence, with the field
%                   binned (true for each binned kind-1 dimension)
%     Rel           relaxation, with fields T1, T2 (µs, Inf: no relaxation)
%                   and on (true/false)
%     productRule   true: multiply the signals of the nuclear subspaces
%     expandFactor  the binning buffers are 2^expandFactor times longer
%                   than the signal, along each binned dimension
%
%   Output:
%     y      signal: scalar (no increments), row vector (1 dimension) or
%            array with one array dimension per kind-1 dimension
%     rates  distinct decay rates (1/µs) of the pathways, one row per group
%            of pathways and one column per kind-1 dimension
%     nOff   number of resonance offsets (largest over the points)
function [y,rates,nOff] = pathwaySignal(E,seq,Rel,productRule,expandFactor)

nOff = 1;
N = seq.N;
nP = seq.nPathways;
nDims = numel(seq.kind1Dims);

% Relaxation: fixed decay factor per pathway and decay rates per dimension;
% pathways with the same rates share a buffer. Only delays with relaxation
% on (Opt.Relaxation) contribute: tRelax, and incSchemeRelax for the
% incremented delays. The decay rate along a dimension is the sum of the
% rates of the incremented intervals, each weighted with how fast it gets
% longer.
coh = seq.co~=0;
rateInt = zeros(nP,N);
if Rel.on
  rateInt(coh) = 1/Rel.T2;
  rateInt(~coh) = 1/Rel.T1;
end
fixedDecay = exp(-rateInt*seq.tRelax.');
pRates = rateInt*seq.incSchemeRelax;
[rates,~,group] = unique(pRates,'rows');
nGroups = size(rates,1);
pref = seq.pathwayprefactor(:).*fixedDecay;

% Stepped dimensions: point by point, with fixed delays at each point
if any(~seq.binned) && ~productRule
  [y,nOff] = steppedSignal(E,seq,Rel,expandFactor);
  return
end

% Binned dimensions: incremented intervals and their incrementation scheme
incIntervals = find(any(seq.incScheme~=0,2)).';
incScheme = seq.incScheme(incIntervals,:);

% Buffers
if nDims>0
  nPointsF = 2^expandFactor*seq.nPointsInc;
  if nDims==1, siz = [1 nPointsF]; else, siz = nPointsF; end
end

nSub = numel(E.Ea);
[t0,tStep] = netCoherenceTimes(seq);
sigma = E.offWidth/(2*sqrt(log(2)));  % 1/e half width of the line (MHz)
if productRule
  % per pathway and subspace, time domain products, per entry
  y = 0;
  stepped = ~any(seq.binned);
  for e = 1:E.nE
    for p = 1:nP
      td_ = 1;
      for s = 1:nSub
        Ed = entrySlice(E,s,e);
        % the electron resonance offset enters once (ideal pulses only)
        if s==1, Ed = addOffsets(Ed,Ed.offCenter); end
        if stepped
          % point by point, with fixed delays at each point
          td1 = zeros([seq.nPointsInc 1 1]);
          for iPt = 1:numel(td1)
            sk = fixDelays(seq,1:nDims,stepsAtPoint(seq.nPointsInc,iPt));
            [~,D] = evaluatePathways(Ed,sk,p,[]);
            td1(iPt) = trace(D);
          end
          if nDims==1, td1 = td1.'; end
          if s==1, td1 = td1*pref(p)*E.w(e); end
          td_ = td_.*td1;
          continue
        end
        [G,D,BL,BR] = evaluatePathways(Ed,seq,p,[]);
        if s==1, G = G*pref(p)*E.w(e); end
        buff = complex(zeros(siz));  % sf_peaks adds to it in place
        sf_peaks(buff,seq.dt,incScheme,seq.idxFreeL(p,incIntervals),seq.idxFreeR(p,incIntervals),...
          Ed.Ea,Ed.Eb,G,D,BL{:},BR{:});
        if nDims==1
          td1 = ifft(buff)*numel(buff);
          td1 = td1(1:seq.nPointsInc);
        else
          td1 = ifft2(buff)*numel(buff);
          td1 = td1(1:seq.nPointsInc(1),1:seq.nPointsInc(2));
        end
        td_ = td_.*td1;
      end
      % exact average over the inhomogeneous line (ideal pulses), with the
      % net electron coherence time Tnet at each point
      Tnet = t0(p);
      if nDims>0, Tnet = Tnet + tStep(p,1)*(0:seq.nPointsInc(1)-1); end
      if nDims>1, Tnet = Tnet.' + tStep(p,2)*(0:seq.nPointsInc(2)-1); end
      td_ = td_.*exp(-(pi*sigma(e)*Tnet).^2);
      td_ = applyDecay(td_,rates(group(p),:),seq);
      y = y + td_;
    end
  end
  return
end

% Without product rule: one subspace
if nDims>0
  % complex buffers, one per group; sf_peaks adds to them in place
  buffers = cell(1,nGroups);
  for iG = 1:nGroups
    buffers{iG} = complex(zeros(siz));
  end
else
  yScalar = zeros(1,nGroups);
end

n = size(E.Ea{1},1);
chunk = max(1,floor(4e6/n^2));
% Average over the inhomogeneous line: exact for ideal pulses at a single
% point (exp(-(pi*sigma*Tnet)^2) per pathway and entry); otherwise over
% the offsets from offsetQuadrature. With ideal pulses, refocused pathways
% do not depend on the offset: they are evaluated once, with weight 1.
exactLineAverage = ~seq.finitePulses && nDims==0;
if exactLineAverage
  offsets = 0; offsetWeights = 1;
else
  [offsets,offsetWeights] = offsetQuadrature(seq,E.offWidth,E.offsetTolerance);
end
nOff = numel(offsets);
offsetFree = ~seq.finitePulses & refocusedPathways(seq);
chunks = entryChunks(E,chunk);
for iChunk = 1:numel(chunks)
  ec = chunks{iChunk};
  Ed0 = entrySlice(E,1,ec);
  for iOff = 1:nOff
    Ed = addOffsets(Ed0,Ed0.offCenter + offsets(iOff)*sigma(ec));
    if seq.finitePulses
      PulseSub = finitePulsePropagators(Ed,seq);
    else
      PulseSub = [];
    end
    memo = struct('keys',{{}},'vals',{{}});
    for p = 1:nP
      if offsetFree(p)
        if iOff>1, continue; end
        wOff = 1;
      else
        wOff = offsetWeights(iOff);
      end
      [G,D,BL,BR,memo] = evaluatePathways(Ed,seq,p,PulseSub,memo);
      wp = pref(p)*wOff*reshape(Ed.w,1,1,[]);
      if exactLineAverage
        wp = wp.*reshape(exp(-(pi*sigma(ec)*t0(p)).^2),1,1,[]);
      end
      if nDims==0
        % no increments: G is the identity, signal = trace(D)
        v = sum(diag3(D).*reshape(wp,1,[]),'all');
        yScalar(group(p)) = yScalar(group(p)) + v;
      else
        G = G.*wp;
        iG = group(p);
        buff = buffers{iG};  % shares its data with buffers{iG}
        % manifolds of the incremented free evolutions of this pathway, and
        % the transfer matrices between them, for one entry at a time
        iL = seq.idxFreeL(p,incIntervals); iR = seq.idxFreeR(p,incIntervals);
        ML = cell(1,numel(BL)); MR = cell(1,numel(BR));
        for k = 1:numel(ec)
          for m = 1:numel(BL)
            ML{m} = BL{m}(:,:,k);
            MR{m} = BR{m}(:,:,k);
          end
          sf_peaks(buff,seq.dt,incScheme,iL,iR,Ed.Ea(:,k),Ed.Eb(:,k),G(:,:,k),D(:,:,k),ML{:},MR{:});
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
  buff = buffers{iG};
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
% chooseBinning  Binned or stepped computation of each kind-1 dimension
%
%   binnedDims = chooseBinning(OptBinning,seq,nStates,nEntries,productRule)
%
%   Along a binned dimension, the peaks of each pathway are binned in the
%   frequency domain (sf_peaks). With K incremented intervals, a pathway
%   has nStates^(2K) peaks per entry, independent of the number of points.
%   Along a stepped dimension, the signal is computed point by point, with
%   fixed delays at each point (steppedSignal); this is exact and costs a
%   pathway evaluation per point.
%
%   With 'auto', the dimensions are stepped if all kind-1 dimensions
%   together have at most 64 points (exact, cheap). Otherwise, of the
%   possible choices (binningProblem), the one with the smallest estimated
%   run time is used (binningCost).
%
%   Input:
%     OptBinning   Opt.Binning: 'auto', true or false for all kind-1
%                  dimensions, or one true/false per dimension of the
%                  experiment (Exp.nPoints)
%     seq          parsed sequence
%     nStates      number of nuclear states (largest subspace)
%     nEntries     number of entries
%     productRule  true: product rule
%
%   Output:
%     binnedDims   true for each binned dimension of the experiment, 1 x
%                  numel(seq.kind)
function binnedDims = chooseBinning(OptBinning,seq,nStates,nEntries,productRule)

nDimsExp = numel(seq.kind);
kind1 = seq.kind1Dims;
binnedDims = false(1,nDimsExp);

if ischar(OptBinning) || isstring(OptBinning)
  if ~strcmp(OptBinning,'auto')
    error('Opt.Binning must be ''auto'', true, false, or one true/false per dimension.');
  end
  if prod(seq.nPointsInc)<=64, return; end  % few points: all stepped (exact)
  % all choices of binned kind-1 dimensions
  best = Inf;
  for c = 0:2^numel(kind1)-1
    binned = logical(bitget(c,1:numel(kind1)));
    if ~isempty(binningProblem(seq,binned,productRule)), continue; end
    cost = binningCost(seq,binned,nStates,nEntries);
    if cost<best
      best = cost;
      binnedDims(:) = false;
      binnedDims(kind1(binned)) = true;
    end
  end
  return
end

% user choice: one value for all kind-1 dimensions, or one per dimension
if ~(islogical(OptBinning) || isnumeric(OptBinning)) || ...
    ~(isscalar(OptBinning) || numel(OptBinning)==nDimsExp)
  error('Opt.Binning must be ''auto'', true, false, or one true/false per dimension (%d).',nDimsExp);
end
if isscalar(OptBinning)
  binnedDims(kind1) = logical(OptBinning);
else
  binnedDims = logical(OptBinning(:).');
  notDelays = find(binnedDims & seq.kind~=1);
  if ~isempty(notDelays)
    error('Opt.Binning: dimension %d does not change delays linearly and cannot be binned.',notDelays(1));
  end
end
msg = binningProblem(seq,binnedDims(kind1),productRule);
if ~isempty(msg)
  error('Opt.Binning: %s.',msg);
end

end

%-------------------------------------------------------------------------------
% binningProblem  Why a choice of binned dimensions is not possible
%
%   msg = binningProblem(seq,binned,productRule)
%
%   Binning is limited to 2 dimensions and to K <= 4 incremented intervals
%   (nStates^(2K) peaks per pathway and entry). With the product rule, all
%   dimensions are binned or none.
%
%   Input:
%     seq          parsed sequence
%     binned       true for each binned kind-1 dimension
%     productRule  true: product rule
%
%   Output:
%     msg          the reason, '' if the choice is possible
function msg = binningProblem(seq,binned,productRule)
msg = '';
K = sum(any(seq.incScheme(:,binned)~=0,2));
if sum(binned)>2
  msg = 'at most 2 dimensions can be binned';
elseif K>4
  msg = sprintf(['binning is limited to 4 incremented intervals, but the binned ' ...
    'dimensions change %d intervals; compute a dimension point by point (false)'],K);
elseif productRule && any(binned) && ~all(binned)
  msg = 'with Opt.ProductRule, all dimensions must be binned or all computed point by point';
end
end

%-------------------------------------------------------------------------------
% binningCost  Estimated run time of a choice of binned dimensions
%
%   cost = binningCost(seq,binned,nStates,nEntries)
%
%   Model of the run time of pathwaySignal per pathway. The signal is
%   evaluated once for each point of the stepped dimensions (once if none
%   is stepped):
%   - evaluating the pathway: cPoint + nEntries*(cEntry0 + cEntry*nStates^3)
%     (matrix products of the transfer matrices)
%   - with binned dimensions, in addition: cBin (Fourier transform) and,
%     per entry, cCall + cPeak*nStates^(2K) for K binned intervals
%   The constants were fitted to saffron run times (2026-11-13; 1D and 2D
%   ESEEM, HYSCORE, Carr-Purcell and 2D-CP with 1 to 4 protons; MATLAB
%   R2024b, Windows). Only their ratios matter for the choice. The model
%   picks the faster method wherever the two differ by more than about 2x.
%
%   Input:
%     seq       parsed sequence
%     binned    true for each binned kind-1 dimension
%     nStates   number of nuclear states
%     nEntries  number of entries
%
%   Output:
%     cost      estimated run time per pathway (s)
function cost = binningCost(seq,binned,nStates,nEntries)
cPoint = 8e-4;    % s, evaluating a pathway at a point
cEntry0 = 7e-6;   % s, per entry
cEntry = 7e-9;    % s, per entry and nStates^3
cBin = 1e-3;      % s, Fourier transform of the binned buffer
cCall = 1.2e-5;   % s, one sf_peaks call (one entry)
cPeak = 8e-9;     % s, one peak
K = sum(any(seq.incScheme(:,binned)~=0,2));
perPoint = cPoint + nEntries*(cEntry0 + cEntry*nStates^3);
if any(binned)
  perPoint = perPoint + cBin + nEntries*(cCall + cPeak*nStates^(2*K));
end
cost = prod(seq.nPointsInc(~binned))*perPoint;
end

%-------------------------------------------------------------------------------
% steppedSignal  Signal along the stepped kind-1 dimensions, point by point
%
%   [y,nOff] = steppedSignal(E,seq,Rel,expandFactor)
%
%   At each point along the stepped dimensions, the incremented delays have
%   fixed lengths, so the sequence is an ordinary one with only the binned
%   dimensions left. Its signal over the binned dimensions is computed by
%   pathwaySignal. Exact along the stepped dimensions.
%
%   Input:
%     E, seq, Rel, expandFactor   see pathwaySignal
%
%   Output:
%     y     signal, with one array dimension per kind-1 dimension, row
%           vector for a single dimension
%     nOff  number of resonance offsets (largest over the points)
function [y,nOff] = steppedSignal(E,seq,Rel,expandFactor)

stepDims = find(~seq.binned);
nPointsStep = seq.nPointsInc(stepDims);
y = zeros([seq.nPointsInc 1 1]);
at = repmat({':'},1,numel(seq.nPointsInc));  % position in y: all binned points
nOff = 1;
for iPt = 1:prod(nPointsStep)
  k = stepsAtPoint(nPointsStep,iPt);
  s = fixDelays(seq,stepDims,k);
  [yPt,~,nOffPt] = pathwaySignal(E,s,Rel,false,expandFactor);
  at(stepDims) = num2cell(k+1);
  y(at{:}) = yPt;
  nOff = max(nOff,nOffPt);
end
if isscalar(seq.nPointsInc), y = y.'; end  % row vector

end

%-------------------------------------------------------------------------------
% fixDelays  Sequence at a point along some kind-1 dimensions
%
%   s = fixDelays(seq,dims,k)
%
%   Lengthens the intervals by k steps along the kind-1 dimensions dims and
%   removes these dimensions: in s, the delays at this point are ordinary
%   delays with fixed lengths.
%
%   Input:
%     seq   parsed sequence
%     dims  kind-1 dimensions (indices into seq.kind1Dims)
%     k     number of steps along each of them, 0-based
%
%   Output:
%     s     sequence with fixed delays along dims
function s = fixDelays(seq,dims,k)
s = seq;
stepTime = k(:).*seq.dt(dims).';  % time (µs) along each dimension
s.t = seq.t + (seq.incScheme(:,dims)*stepTime).';
s.tRelax = seq.tRelax + (seq.incSchemeRelax(:,dims)*stepTime).';
keep = setdiff(1:numel(seq.kind1Dims),dims);
s.kind1Dims = seq.kind1Dims(keep);
s.incScheme = seq.incScheme(:,keep);
s.incSchemeRelax = seq.incSchemeRelax(:,keep);
s.dt = seq.dt(keep);
s.nPointsInc = seq.nPointsInc(keep);
if isfield(seq,'binned'), s.binned = seq.binned(keep); end
end

%-------------------------------------------------------------------------------
% stepsAtPoint  Number of steps along each dimension at a point of a grid
%
%   k = stepsAtPoint(nPoints,iPt)
%
%   Input:
%     nPoints  number of points along each dimension
%     iPt      linear index of the point (1-based, first dimension fastest)
%
%   Output:
%     k        number of steps along each dimension, 0-based, row vector
function k = stepsAtPoint(nPoints,iPt)
k = cell(1,numel(nPoints));
if isempty(k), k = zeros(1,0); return; end
[k{:}] = ind2sub([nPoints 1],iPt);
k = [k{:}]-1;
end

%-------------------------------------------------------------------------------
% applyDecay  Apply relaxation decays along the incremented dimensions
%
%   td = applyDecay(td,rate,seq)
%
%   Input:
%     td    time-domain signal, row vector (1 dimension) or matrix (2 dimensions)
%     rate  decay rate (1/µs) of each dimension
%     seq   parsed sequence; uses dt, nPointsInc
%
%   Output:
%     td    signal with the decays applied
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
% entrySlice  Extract the data of some entries for one subspace
%
%   Ed = entrySlice(E,s,ec)
%
%   Input:
%     E   entries, see the header of sf_engine_manifold
%     s   subspace index
%     ec  entry indices
%
%   Output:
%     Ed  data of the entries, in page form (n x 1 x K, n x n x K), see the
%         header of sf_engine_manifold
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
Ed.offCenter = E.offCenter(ec);
Ed.offWidth = E.offWidth(ec);
if isempty(E.moment)
  Ed.moment = [];
else
  Ed.moment = E.moment(ec);
end
Ed.n = n;
Ed.K = K;
% For field-reversed partners (E.rev, appended after the original entries):
% position of the original entry in this slice, 0 if not present
Ed.original = zeros(1,K);
rev = E.rev(ec);
if any(rev)
  [found,loc] = ismember(ec(rev)-E.nE/2,ec);
  idx = find(rev);
  Ed.original(idx(found)) = loc(found);
end
end

%-------------------------------------------------------------------------------
% entryChunks  Divide the entries into chunks
%
%   chunks = entryChunks(E,chunkSize)
%
%   Field-reversed partners go into the same chunk as their original
%   entries.
%
%   Input:
%     E          entries, see the header of sf_engine_manifold
%     chunkSize  maximum number of entries per chunk
%
%   Output:
%     chunks     cell array of vectors of entry indices
function chunks = entryChunks(E,chunkSize)
if any(E.rev)
  nOriginal = E.nE/2;
  chunkSize = max(1,floor(chunkSize/2));
  starts = 1:chunkSize:nOriginal;
  chunks = cell(1,numel(starts));
  for c = 1:numel(starts)
    ec = starts(c):min(starts(c)+chunkSize-1,nOriginal);
    chunks{c} = [ec ec+nOriginal];
  end
else
  starts = 1:chunkSize:E.nE;
  chunks = cell(1,numel(starts));
  for c = 1:numel(starts)
    chunks{c} = starts(c):min(starts(c)+chunkSize-1,E.nE);
  end
end
end

%-------------------------------------------------------------------------------
% addOffsets  Shift entry energies by their electron resonance offsets
%
%   Ed = addOffsets(Ed,offsets)
%
%   Gives the entry data in the frame rotating at Exp.mwFreq. Free evolution
%   and pulses both use the shifted energies.
%
%   Input:
%     Ed       entry data, from entrySlice
%     offsets  electron resonance offset of each entry (MHz, including the
%              offset within the inhomogeneous line), 1 x K
%
%   Output:
%     Ed       entry data, with the energies of the two manifolds shifted
%              by +offsets/2 (a) and -offsets/2 (b)
function Ed = addOffsets(Ed,offsets)
Ed.Ea = Ed.Ea + offsets/2;
Ed.Eb = Ed.Eb - offsets/2;
Ed.Ea3 = reshape(Ed.Ea,Ed.n,1,Ed.K);
Ed.Eb3 = reshape(Ed.Eb,Ed.n,1,Ed.K);
end

%-------------------------------------------------------------------------------
% refocusedPathways  Find the pathways that refocus at the detection point
%
%   refocused = refocusedPathways(seq)
%
%   A pathway refocuses at the detection point, for every point of the
%   incremented delays, if it has zero net electron coherence time over the
%   free evolution (also during RF pulses) and in each incremented
%   dimension. With ideal pulses, its signal does not depend on the
%   resonance offset.
%
%   Input:
%     seq        parsed sequence; uses t and the fields used by
%                netCoherenceTimes
%
%   Output:
%     refocused  true for each refocused pathway, nPathways x 1
function refocused = refocusedPathways(seq)
[t0,tStep] = netCoherenceTimes(seq);
refocused = abs(t0)<=1e-12*max([1 seq.t]) & all(tStep==0,2);
end

%-------------------------------------------------------------------------------
% netCoherenceTimes  Net electron coherence time of each pathway
%
%   [t0,tStep] = netCoherenceTimes(seq)
%
%   The net electron coherence time of a pathway is the sum of the free
%   evolution times (also during RF pulses), each multiplied by the
%   electron coherence order of its interval. It changes linearly with the
%   increments of the incremented delays.
%
%   Input:
%     seq    parsed sequence; uses co, t, isRF, Pulse, incScheme, dt, nPointsInc
%
%   Output:
%     t0     net coherence time (µs) at the first point, nPathways x 1
%     tStep  change (µs) per increment, nPathways x nDims
function [t0,tStep] = netCoherenceTimes(seq)
tFree = seq.t + seq.isRF.*[seq.Pulse.tp];
t0 = seq.co*tFree(:);
nDims = numel(seq.nPointsInc);
tStep = zeros(numel(t0),nDims);
for d = 1:nDims
  tStep(:,d) = seq.co*seq.incScheme(:,d)*seq.dt(d);
end
end

%-------------------------------------------------------------------------------
% offsetQuadrature  Offsets and weights for the average over the EPR line
%
%   [x,w] = offsetQuadrature(seq,offWidth,tol)
%
%   Gives the resonance offsets and weights for averaging the signal over
%   the inhomogeneous EPR line, with as few offsets as the accuracy tol
%   allows for the timing of the sequence. If the signal
%   does not depend on the offset (no inhomogeneous line, or ideal pulses
%   with all pathways refocused at the detection point), it gives a single
%   offset.
%
%   Method: the offsets of each entry have a Gaussian distribution,
%   offCenter + x*sigma with weight exp(-x^2), where sigma =
%   offWidth/(2*sqrt(log(2))) is the 1/e half width. The trapezoid rule
%   samples x over |x| <= X with step h. X and h are chosen so that each of
%   the two errors is of the order of tol times the signal amplitude:
%   - Truncation: the weight beyond X is erfc(X), so X = erfcinv(tol).
%   - Aliasing: the step h*sigma (MHz) adds copies of the signal of each
%     pathway at net electron coherence times shifted by 1/(h*sigma). The
%     signal extends over |t| <= Teff, the longest net coherence time of
%     the free evolution plus the total length of the finite pulses (whose
%     propagators vary with the offset on the scale of 1/tp). The decay
%     exp(-(pi*sigma*t)^2) of the line widens this to Teff + c/(pi*sigma),
%     with c = sqrt(log(1/tol)). So h <= 1/(sigma*Teff + c/pi), for the
%     largest sigma.
%
%   Input:
%     seq       parsed sequence; uses finitePulses, isRF, Pulse, nPointsInc
%               and the fields used by netCoherenceTimes
%     offWidth  FWHM (MHz) of the offset distribution of each entry, 1 x nE
%     tol       error target (Opt.OffsetTolerance)
%
%   Output:
%     x    offsets in units of sigma, 1 x nOff
%     w    weights, 1 x nOff, sum 1
function [x,w] = offsetQuadrature(seq,offWidth,tol)
sigma = max([offWidth 0])/(2*sqrt(log(2)));  % largest 1/e half width (MHz)
offsetDependent = seq.finitePulses | ~refocusedPathways(seq);
if sigma==0 || ~any(offsetDependent)
  x = 0; w = 1;
  return
end

% Net electron coherence time of each pathway: from t0 at the first point,
% it changes linearly by dT over the range of each incremented dimension,
% so its extremes are tMin and tMax
[t0,tStep] = netCoherenceTimes(seq);
t0 = t0(offsetDependent);
dT = tStep(offsetDependent,:).*(seq.nPointsInc(:).'-1);
tMin = t0 + sum(min(dT,0),2);
tMax = t0 + sum(max(dT,0),2);
tp = [seq.Pulse.tp];
Teff = max(abs([tMin; tMax])) + sum(tp(~seq.isRF & ~[seq.Pulse.ideal]));

% Trapezoid rule
X = erfcinv(tol);
c = sqrt(log(1/tol));
nHalf = ceil(X*(sigma*Teff + c/pi));  % from h <= 1/(sigma*Teff + c/pi)
x = linspace(-X,X,2*nHalf+1);
w = exp(-x.^2);
w([1 end]) = w([1 end])/2;
w = w/sum(w);  % offset-independent signals exact (error erfc(X) otherwise)
end

%-------------------------------------------------------------------------------
% finitePulsePropagators  Sub-propagators of the finite microwave pulses
%
%   PulseSub = finitePulsePropagators(Ed,seq)
%
%   A pulse with frequency offset f (Pulse.Frequency) from Exp.mwFreq is
%   static in the frame rotating at mwFreq + f, where the Hamiltonian is
%   Hp = H0 - f*Sz + pulse term, with the pulse phase referenced to the
%   start of the pulse (as in pulse.m). The propagator is computed in that
%   frame and transformed back to the frame rotating at mwFreq, in which the
%   free evolution and the detection are computed:
%   U = exp(-2i*pi*f*tp*Sz)*exp(-2i*pi*Hp*tp). Without pulse term (Flip = 0),
%   this is the free evolution exp(-2i*pi*H0*tp).
%
%   Input:
%     Ed        entry data, with the energies including the resonance
%               offsets (see addOffsets)
%     seq       parsed sequence; uses N, isRF, Pulse
%
%   Output:
%     PulseSub  sub-propagators for all entries, n x n x K, in a cell array
%               with one row per interval and one column per block (1 = aa,
%               2 = bb, 3 = ab, 4 = ba); empty for RF and ideal pulses
function PulseSub = finitePulsePropagators(Ed,seq)
N = seq.N;
n = Ed.n; K = Ed.K;
PulseSub = cell(N,4);
idxa = 1:n; idxb = n+1:2*n;
Szdiag = [ones(n,1); -ones(n,1)]/2;
for i = 1:N
  if seq.isRF(i) || seq.Pulse(i).ideal, continue; end
  tp = seq.Pulse(i).tp;
  % nutation frequency of each transition (MHz): pulse amplitude times the
  % transition moment |<a|S+|b>| (1 for S = 1/2 with all nuclei)
  nu1 = seq.Pulse(i).Flip/(2*pi*tp)*ones(1,K);
  if ~isempty(Ed.moment), nu1 = nu1.*Ed.moment; end
  phaseFactor = exp(-1i*seq.Pulse(i).Phase);
  f = seq.Pulse(i).Frequency;  % MHz, offset from Exp.mwFreq
  % propagators in the frame of the pulse
  U = zeros(2*n,2*n,K);
  for k = 1:K
    j = Ed.original(k);
    if j>0 && isequal(Ed.Ea(:,k),Ed.Ea(:,j)) && isequal(Ed.Eb(:,k),Ed.Eb(:,j))
      % Field-reversed partner: its pulse Hamiltonian is conj(Hp) = Hp.' of
      % the original entry, so its propagator is the transpose
      U(:,:,k) = U(:,:,j).';
      continue
    end
    Hp = [diag(Ed.Ea(:,k)), Ed.M(:,:,k)*nu1(k)/2; Ed.Mt(:,:,k)*nu1(k)/2, diag(Ed.Eb(:,k))];
    Hp = Hp - diag(f*Szdiag);
    [V,lam] = eig((Hp+Hp')/2,'vector');
    U(:,:,k) = (V.*exp(-2i*pi*tp*lam.'))*V';
  end
  % back to the frame rotating at Exp.mwFreq: exp(-2i*pi*f*tp*Sz)*U, for
  % all entries (after the transposes of the field-reversed partners, whose
  % pulse-frame propagators are the transposes)
  if f~=0
    U = exp(-2i*pi*f*tp*Szdiag).*U;
  end
  PulseSub{i,1} = U(idxa,idxa,:);
  PulseSub{i,2} = U(idxb,idxb,:);
  PulseSub{i,3} = U(idxa,idxb,:)*phaseFactor;
  PulseSub{i,4} = U(idxb,idxa,:)*conj(phaseFactor);
end
end

%-------------------------------------------------------------------------------
% evaluatePathways  Evaluate one coherence transfer pathway for all entries
%
%   [G,D,BL,BR,memo] = evaluatePathways(Ed,seq,p,PulseSub)
%   [G,D,BL,BR,memo] = evaluatePathways(Ed,seq,p,PulseSub,memo)
%
%   Operands are diagonal (n x 1 x K) or full (n x n x K). Products of
%   prefixes are shared between pathways (memo).
%
%   Input:
%     Ed        entry data, from entrySlice
%     seq       parsed sequence
%     p         pathway index
%     PulseSub  finite pulse sub-propagators, from finitePulsePropagators;
%               [] for ideal pulses
%     memo      cache of the products of earlier pathways, optional, see
%               memoGet
%
%   Output:
%     G         state before the first incremented delay, n x n x K; empty
%               without incremented delays
%     D         detection operator S+, propagated back from the end,
%               n x n x K (without incremented delays, trace(D) is the
%               signal)
%     BL, BR    left and right transfer blocks between consecutive
%               incremented delays, cell arrays of n x n x K arrays
%     memo      updated cache
function [G,D,BL,BR,memo] = evaluatePathways(Ed,seq,p,PulseSub,memo)
if nargin<5 || isempty(memo), memo = struct('keys',{{}},'vals',{{}}); end
N = seq.N;
incIntervals = find(any(seq.incScheme~=0,2)).';
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

% memoGet  Look up a value in the cache of evaluatePathways
%
%   [hit,val] = memoGet(memo,key)
%
%   Input:
%     memo  cache, structure with cell array fields keys and vals
%     key   key string, identifying a product of propagator blocks
%
%   Output:
%     hit   true if the key is in the cache
%     val   cached value, empty if the key is not in the cache
function [hit,val] = memoGet(memo,key)
k = find(strcmp(memo.keys,key),1);
hit = ~isempty(k);
if hit, val = memo.vals{k}; else, val = []; end
end

% memoSet  Add a value to the cache of evaluatePathways
%
%   memo = memoSet(memo,key,val)
%
%   Input:
%     memo  cache, see memoGet
%     key   key string
%     val   value to store
%
%   Output:
%     memo  updated cache
function memo = memoSet(memo,key,val)
memo.keys{end+1} = key;
memo.vals{end+1} = val;
end

% stepBlock  Propagator block of one interval
%
%   S = stepBlock(Ed,seq,PulseSub,i,c)
%
%   An interval is a pulse, followed by free evolution over the initial
%   delay. The free evolution is in the manifold of the row index.
%
%   Input:
%     Ed        entry data, from entrySlice
%     seq       parsed sequence; uses t, isRF, Pulse
%     PulseSub  finite pulse sub-propagators, from finitePulsePropagators
%     i         interval index
%     c         block index: 1 = aa, 2 = bb, 3 = ab, 4 = ba
%
%   Output:
%     S         propagator block, diagonal (n x 1 x K) or full (n x n x K)
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

% mulLeft  Compute P*L for stacks of matrices
%
%   X = mulLeft(P,L)
%
%   Input:
%     P  stack of matrices, diagonal (n x 1 x K) or full (n x n x K)
%     L  stack of matrices, empty (identity), diagonal or full
%
%   Output:
%     X  product P*L, diagonal or full
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

% mulRight  Compute R*P' for stacks of matrices
%
%   X = mulRight(R,P)
%
%   Input:
%     R  stack of matrices, empty (identity), diagonal (n x 1 x K) or full
%        (n x n x K)
%     P  stack of matrices, diagonal or full
%
%   Output:
%     X  product R*P', diagonal or full
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

% mulFull  Compute A*B for stacks of matrices, with a full result
%
%   X = mulFull(A,B,n)
%
%   Input:
%     A, B  stacks of matrices, empty (identity), diagonal (n x 1 x K) or
%           full (n x n x K)
%     n     matrix dimension
%
%   Output:
%     X     product A*B, full; empty if both A and B are empty
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

% toFull  Convert a stack of matrices to full form
%
%   X = toFull(A,n)
%
%   Input:
%     A  stack of matrices, empty (identity), diagonal (n x 1 x K) or full
%        (n x n x K)
%     n  matrix dimension
%
%   Output:
%     X  full matrices, n x n x K
function X = toFull(A,n)
if isempty(A)
  X = eye(n);
elseif size(A,2)==1 && n>1
  X = eye(n).*permute(A,[2 1 3]);
else
  X = A;
end
end

% diag3  Diagonals of a stack of matrices
%
%   d = diag3(X)
%
%   Input:
%     X  stack of matrices, n x n x K
%
%   Output:
%     d  diagonals of the matrices, n x K
function d = diag3(X)
n = size(X,1);
K = size(X,3);
X = reshape(X,n*n,K);
d = X(1:n+1:end,:);
end


%===============================================================================
% fieldReversalShortcut  Shortcut for the field-reversed partners
%
%   shortcut = fieldReversalShortcut(seq)
%
%   Determines how to obtain the field-reversed partners of a partial-sphere
%   grid. Reversing the field keeps the energies and turns M into conj(M),
%   so each propagator becomes its transpose (conj(H) = H.' for Hermitian
%   H). The reversed signal is therefore that of the same orientation for
%   the time-reversed sequence, V(-n,seq) = V(n,seq'). If seq' is equivalent
%   to seq (e.g. 2p, 3p, 4p ESEEM with ideal pulses, square HYSCORE), the
%   reversed signal follows from the signal itself. The check evaluates a
%   random small system with M and conj(M) at a few points, without binning.
%
%   Input:
%     seq       parsed sequence, from sf_parsesequence
%
%   Output:
%     shortcut  "conj": the reversed signal is the complex conjugate of the
%               signal; "conjtranspose": it is the conjugate transpose;
%               "compute": the partners have to be computed
function shortcut = fieldReversalShortcut(seq)

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

% pointSignal  Signal at one point of the incremented dimensions
%
%   v = pointSignal(Ed,seq,k)
%
%   Evaluates the signal without binning.
%
%   Input:
%     Ed   entry data of a single entry
%     seq  parsed sequence
%     k    number of increments along each incremented dimension, 0-based
%
%   Output:
%     v    signal at this point
function v = pointSignal(Ed,seq,k)
s = fixDelays(seq,1:numel(seq.kind1Dims),k);
v = 0;
for p = 1:seq.nPathways
  [~,D] = evaluatePathways(Ed,s,p,[]);
  v = v + seq.pathwayprefactor(p)*trace(D);
end
end


%===============================================================================
% endorSweep  ENDOR signal as a function of the RF frequency
%
%   [y,nOff,nDirections] = endorSweep(E,seq,NucHams,nSubSpaces,rf)
%
%   Computes the echo over the RF frequency (kind-2 dimension), for a single
%   rectangular RF pulse. The RF pulse transfers population between the two
%   levels u,v of each NMR transition with probability
%   P(Delta) = (nu2/Omega)^2*sin(pi*Omega*tp)^2, Omega = sqrt(nu2^2+Delta^2),
%   at the detuning Delta = nu_uv - nu_RF from the absolute RF frequency
%   nu_RF (carrier plus pulse frequency, seq.rfAxis). Terms that depend on
%   the RF phase average out (random RF phase), leaving the population
%   transfer.
%
%   RF field direction: the linearly polarized RF field B2 along a
%   direction n in the plane perpendicular to the static field gives each
%   transition the nutation frequency nu2 = 2*B2*|<u|O(n)|v>| (rfNutation).
%   Crystals have n = yL, perpendicular to the microwave field (along xL).
%   In powders, the signal is averaged over n = cos(chi)*xL + sin(chi)*yL,
%   chi in [0,pi) (see sf_rfdirectionaverage). The average is over the
%   signal, since P is nonlinear in nu2.
%
%   Input:
%     E           entries, see the header of sf_engine_manifold
%     seq         parsed sequence, from sf_parsesequence
%     NucHams     nuclear sub-Hamiltonians, from nuclearHamiltonians
%     nSubSpaces  number of nuclear subspaces (must be 1)
%     rf          RF pulse, with fields
%       B2          RF field amplitude (mT)
%       tp          pulse length (µs)
%       averageDirections  true to average over the RF field direction
%                   (powders), false for the RF field along yL only
%                   (crystals)
%
%   Output:
%     y            signal over the RF axis (seq.rfAxis), row vector
%     nOff         number of resonance offsets
%     nDirections  number of RF field directions in the average
function [y,nOff,nDirections] = endorSweep(E,seq,NucHams,nSubSpaces,rf)

if nSubSpaces>1
  error('ENDOR is not supported with the product rule.');
end
iRF = find(seq.isRF);
if iRF==1
  error('The RF pulse cannot be the first pulse.');
end
% The population-transfer model needs the electron spin in a population
% state during the RF pulse (blocks 1 = aa, 2 = bb, not 3 = ab, 4 = ba)
if any(seq.ctp(:,iRF)>2)
  error(['The manifold engine simulates ENDOR only for pathways with electron coherence order 0 ' ...
    'during the RF pulse. Select them with Exp.Pathways (e.g. ''+00-'' for Mims ENDOR), ' ...
    'or use Opt.Method = ''full''.']);
end

% Sequence split at the RF pulse: intervals 1..iRF-1 (state G before the RF
% pulse) and iRF..N (detection D after it). G = state at the end of
% interval iRF-1 (evaluatePathways splits the pathway after each
% incremented interval)
seqG = seq;
seqG.incScheme = zeros(seq.N,1);
seqG.incScheme(iRF-1) = 1;

n = size(E.Ea{1},1);
chunk = max(1,floor(4e6/n^2));
% Average over the inhomogeneous line: exact for ideal pulses
% (exp(-(pi*sigma*Tnet)^2) per pathway and entry), otherwise over the
% offsets from offsetQuadrature
sigma = E.offWidth/(2*sqrt(log(2)));  % 1/e half width of the line (MHz)
if seq.finitePulses
  [offsets,offsetWeights] = offsetQuadrature(seq,E.offWidth,E.offsetTolerance);
else
  offsets = 0; offsetWeights = 1;
  Tnet = netCoherenceTimes(seq);
end
nOff = numel(offsets);
[idxU,idxV] = find(triu(true(n),1));  % levels u<v of each NMR transition
nTransitions = numel(idxU);

% Signal without RF pulse (yBaseline), and population transfer amplitudes,
% NMR frequencies and RF matrix elements (for fields along xL and yL) of
% the NMR transitions of all entries (amp, nuNMR, ox, oy, one cell per
% chunk of entries)
yBaseline = 0;
chunks = entryChunks(E,chunk);
nChunks = numel(chunks);
[amp,nuNMR,ox,oy] = deal(cell(1,nChunks));
for iChunk = 1:nChunks
  ec = chunks{iChunk};
  Ed0 = entrySlice(E,1,ec);
  K = numel(ec);
  % RF matrix elements (MHz/mT) of the NMR transitions, alpha manifold in
  % the first nTransitions rows, beta in the rest
  [oxc,oyc] = rfNutation(E,ec,NucHams,idxU,idxV);
  % population transfer amplitudes of the NMR transitions in the manifold
  % of the RF pulse, summed over offsets and pathways, rows as in oxc
  transferAmp = zeros(2*nTransitions,K);
  for iOff = 1:nOff
    Ed = addOffsets(Ed0,Ed0.offCenter + offsets(iOff)*sigma(ec));
    if seq.finitePulses
      PulseSub = finitePulsePropagators(Ed,seq);
    else
      PulseSub = [];
    end
    for p = 1:seq.nPathways
      [G,D] = evaluatePathways(Ed,seqG,p,PulseSub);
      wp = seq.pathwayprefactor(p)*offsetWeights(iOff)*reshape(Ed.w,1,1,[]);
      if ~seq.finitePulses
        wp = wp.*reshape(exp(-(pi*sigma(ec)*Tnet(p)).^2),1,1,[]);
      end
      G = G.*wp;
      % baseline: trace(D*G)
      yBaseline = yBaseline + sum(G.*permute(D,[2 1 3]),'all');
      % population transfer
      popG = diag3(G); popD = diag3(D);
      rows = (seq.ctp(p,iRF)-1)*nTransitions + (1:nTransitions);  % manifold of the RF pulse
      transferAmp(rows,:) = transferAmp(rows,:) + ...
        (popG(idxV,:)-popG(idxU,:)).*(popD(idxU,:)-popD(idxV,:));
    end
  end
  nuNMRc = abs([Ed0.Ea(idxV,:)-Ed0.Ea(idxU,:); Ed0.Eb(idxV,:)-Ed0.Eb(idxU,:)]);
  keep = transferAmp~=0;
  amp{iChunk} = transferAmp(keep);
  nuNMR{iChunk} = nuNMRc(keep);
  ox{iChunk} = oxc(keep);
  oy{iChunk} = oyc(keep);
end
amp = vertcat(amp{:}); nuNMR = vertcat(nuNMR{:});
ox = vertcat(ox{:}); oy = vertcat(oy{:});
maxTransferAmp = max([0; abs(amp)]);  % largest population transfer amplitude

% RF response: depends only on the NMR transition frequencies and nutation
% frequencies, not on the offset or the pathway
keep = ox~=0 | oy~=0;  % no RF matrix element: no transfer
amp = amp(keep); nuNMR = nuNMR(keep); ox = ox(keep); oy = oy(keep);
rfFreqs = seq.rfAxis;  % absolute RF frequencies (MHz)
nu2 = @(chi) 2*rf.B2*abs(cos(chi)*ox + sin(chi)*oy);  % RF field along cos(chi)*xL + sin(chi)*yL
if ~rf.averageDirections
  nDirections = 1;
  y = yBaseline + rfResponse(amp,nuNMR,2*rf.B2*abs(oy),rfFreqs,rf.tp);  % RF field along yL
else
  signal = @(chi) yBaseline + rfResponse(amp,nuNMR,nu2(chi),rfFreqs,rf.tp);
  [y,nDirections,converged] = sf_rfdirectionaverage(signal,pi);
  if ~converged
    warning(['The average over the RF field direction did not converge with %d directions. ' ...
      'The result can be inaccurate.'],nDirections);
  end
end

% No population transfer by the RF pulse at any RF frequency, although there
% is an echo: the ENDOR spectrum is flat. This happens for symmetric
% excitation in Davies ENDOR (inversion pulse centered between the two
% hyperfine lines of a spin-1/2 nucleus, with no inhomogeneous line width),
% and at the blind spots of Mims ENDOR.
if yBaseline~=0 && maxTransferAmp<=1e-12*abs(yBaseline)
  warning('EasySpin:saffron:ENDORFlat',...
    ['The RF pulse does not change the echo, so the ENDOR spectrum is flat. In Davies ENDOR, ' ...
     'this happens if the inversion pulse inverts the hyperfine lines of the nucleus equally, ' ...
     'e.g. without Exp.mwFreq and Sys.lwAnisoG, where the pulses are centered on the EPR line. ' ...
     'Give an inhomogeneous line width in Sys.lwAnisoG or an off-center Exp.mwFreq. ' ...
     'In Mims ENDOR, it happens at the blind spots (tau = k/A).']);
end
end

%-------------------------------------------------------------------------------
% rfNutation  RF matrix elements of the NMR transitions
%
%   [ox,oy] = rfNutation(E,ec,NucHams,idxU,idxV)
%
%   Computes the matrix elements <u|O(n)|v> of the RF operator O(n)
%   (MHz/mT) for an RF field along n = xL and n = yL, in the eigenbases of
%   the nuclear Hamiltonians of the two manifolds. O(n) is the nuclear
%   Zeeman interaction along n (with gnscale and chemical shielding) plus
%   the hyperfine enhancement dS(n)*A*I, from the tilt dS(n) of the
%   electron spin by the RF field (E.dSa, E.dSb, zero without
%   enhancement). O is linear in n, so a field B2 along
%   cos(chi)*xL + sin(chi)*yL gives the nutation frequency
%   nu2 = 2*B2*|cos(chi)*ox + sin(chi)*oy|.
%
%   Input:
%     E           entries, see the header of sf_engine_manifold
%     ec          entry indices
%     NucHams     nuclear sub-Hamiltonians, from nuclearHamiltonians
%     idxU, idxV  levels u and v of each NMR transition
%
%   Output:
%     ox, oy      matrix elements (MHz/mT), 2*nTransitions x K: alpha
%                 manifold in the first nTransitions rows, beta in the rest
function [ox,oy] = rfNutation(E,ec,NucHams,idxU,idxV)
K = numel(ec);
nTransitions = numel(idxU);
ox = zeros(2*nTransitions,K);
oy = zeros(2*nTransitions,K);
N = NucHams(1);
idx = sub2ind(size(N.Hnzz),idxU,idxV);
for k = 1:K
  e = ec(k);
  xL = E.xL(:,e); yL = E.yL(:,e);
  for iManifold = 1:2
    if iManifold==1
      V = E.Va{1}(:,:,e); dS = E.dSa(:,e);
    else
      V = E.Vb{1}(:,:,e); dS = E.dSb(:,e);
    end
    Ox = xL(1)*N.Hnzx + xL(2)*N.Hnzy + xL(3)*N.Hnzz + ...
      dS(1)*N.Hhfx + dS(2)*N.Hhfy + dS(3)*N.Hhfz;
    Oy = yL(1)*N.Hnzx + yL(2)*N.Hnzy + yL(3)*N.Hnzz + ...
      dS(4)*N.Hhfx + dS(5)*N.Hhfy + dS(6)*N.Hhfz;
    OxEig = V'*Ox*V;
    OyEig = V'*Oy*V;
    rows = (iManifold-1)*nTransitions+(1:nTransitions);
    ox(rows,k) = OxEig(idx);
    oy(rows,k) = OyEig(idx);
  end
end
end

%-------------------------------------------------------------------------------
% rfResponse  Signal change by population transfer of a rectangular RF pulse
%
%   y = rfResponse(amp,nuNMR,nu2,rfFreqs,tp)
%
%   Input:
%     amp      population transfer amplitude of each NMR transition, column
%     nuNMR    NMR frequency (MHz) of each transition, column
%     nu2      RF nutation frequency (MHz) of each transition, column
%     rfFreqs  absolute RF frequencies (MHz), row
%     tp       RF pulse length (µs)
%
%   Output:
%     y        sum over the transitions of amp*P, with the transfer
%              probability P = (nu2/Omega)^2*sin(pi*Omega*tp)^2 at the
%              detuning Delta = nuNMR - rfFreq, Omega = sqrt(nu2^2+Delta^2),
%              row over rfFreqs
function y = rfResponse(amp,nuNMR,nu2,rfFreqs,tp)
y = zeros(1,numel(rfFreqs));
for iFreq = 1:numel(rfFreqs)
  Delta = nuNMR - rfFreqs(iFreq);
  Omega = sqrt(nu2.^2 + Delta.^2);
  P = (nu2./Omega).^2.*sin(pi*Omega*tp).^2;
  P(Omega==0) = 0;
  y(iFreq) = sum(amp.*P);
end
end
