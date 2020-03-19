function [T,F,System,Symmetry,isFieldDep] = magint(System,SpinOps,CenterField,...
  includeNuclearZeeman,includeNuclearQuadrupole,explicitFieldSweep)
% generate the irreducible spherical tensor components and, if using the
% general Liouvillian method, their corresponding operators as well

if isempty(SpinOps)
  % no need to use spin operators for fast code
  generalLiouvillian = false;
else
  generalLiouvillian = true;
end

% Transformation from molecular frame to diffusion frame
% (DiffFrame contains Euler angles for mol->Diff transformation)
R_M2Diff = erot(System.DiffFrame);

% Count the number of interaction terms in the spin Hamiltonian
% -------------------------------------------------------------------------
nElSpins = System.nElectrons;
System.I = nucspin(System.Nucs);
nNucSpins = System.nNuclei;

nElZeeman = nElSpins;
nHyperfine = nElSpins*nNucSpins;
nZFS = sum(System.S>1/2);
nElPairs = nElSpins*(nElSpins-1)/2;
nInteractions = nElZeeman + nHyperfine + nZFS + nElPairs;

if includeNuclearZeeman
  nInteractions = nInteractions + nNucSpins;
end

if includeNuclearQuadrupole
  nInteractions = nInteractions + nNucSpins;
end

T0 = cell(nInteractions,1);
T1 = cell(nInteractions,3);
T2 = cell(nInteractions,5);
F0 = zeros(nInteractions,1);
F1 = zeros(nInteractions,3);
F2 = zeros(nInteractions,5);
isFieldDep = zeros(nInteractions,1);

iInt = 1;

if explicitFieldSweep
  B0 = {0 0 1}; % per Tesla
else
  B0 = {0 0 CenterField/1e3}; % mT -> T
end

% Electron Zeeman interaction terms (muB*B*g*S/h)
%--------------------------------------------------------------------------
for iElSpin = 1:nElSpins
  if ~System.fullg
    g = diag(System.g(iElSpin,:));
  else
    g = System.g((iElSpin-1)*3+(1:3),:);
  end
  if isfield(System,'gFrame')
    R_g2M = erot(System.gFrame(iElSpin,:)).';
    g = R_g2M*g*R_g2M.';  % g frame -> molecular frame
  end
  if ~generalLiouvillian
    % this code wasn't originally in magint, so does the general
    % Liouvillian scheme process a diffusion tilt elsewhere?
    g = R_M2Diff*g*R_M2Diff.';  % molecular frame -> diffusion frame
  end
  if generalLiouvillian
    [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(B0,SpinOps(iElSpin,:));
  end
  [F0(iInt),F1(iInt,:),F2(iInt,:)] = tensor_cart2sph(g*bmagn/planck); % Hz
  isFieldDep(iInt) = true;
  if ~generalLiouvillian
    % Generate custom fields for fast method
    if any(abs(F1(iInt,:)/F0(iInt))>1e-9)
      error('g tensor must be symmetric for this method.');
    end
    % Set parameters for chili_liouvmatrix*
    System.g_axial(iElSpin) = F2(iInt,1)==0;
    System.EZ0(iElSpin) = F0(iInt)*B0{3}*2*pi; % Hz -> angular frequency
    System.EZ2(:,iElSpin) = F2(iInt,:).'*B0{3}*2*pi; % Hz -> angular frequency
  end
  iInt = iInt + 1;
end

% Hyperfine interaction terms (S*A*I)
%--------------------------------------------------------------------------
for iElSpin = 1:nElSpins
  eidx = 3*(iElSpin-1)+(1:3);
  for iNucSpin = 1:nNucSpins
    if System.fullA
      A_ = System.A(3*(iNucSpin-1)+(1:3),eidx)*1e6; % MHz -> Hz
    else
      A_ = diag(System.A(iNucSpin,eidx))*1e6; % MHz -> Hz
    end
    if any(System.AFrame(iNucSpin,eidx))
      R_M2A = erot(System.AFrame(iNucSpin,eidx)); % mol frame -> A frame
      R_A2M = R_M2A.'; % A frame -> mol frame
      A_ = R_A2M*A_*R_A2M.';
    end
    if ~generalLiouvillian
      % this code wasn't originally in magint, so does the general
      % Liouvillian scheme process a diffusion tilt elsewhere?
      A_ = R_M2Diff*A_*R_M2Diff.';  % molecular frame -> diffusion frame
    end
    if generalLiouvillian
      S_ = SpinOps(iElSpin,:);
      I_ = SpinOps(nElSpins+iNucSpin,:);
      [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(S_,I_);
    end
    [F0(iInt),F1(iInt,:),F2(iInt,:)] = tensor_cart2sph(A_);
    isFieldDep(iInt) = false;
    if ~generalLiouvillian
      % Generate custom fields for fast method
      if (any(abs(F1(iInt,:))>1e-6))
        error('Hyperfine tensors must be symmetric for this method.');
      end
      % Set parameters for chili_liouvmatrix*
      System.A_axial(iNucSpin,iElSpin) = F2(iInt,1)==0;
      System.HF0(iNucSpin,iElSpin) = F0(iInt)*2*pi; % Hz -> angular frequency
      System.HF2(:,iNucSpin,iElSpin) = F2(iInt,:).'*2*pi; % Hz -> angular frequency
    end
    iInt = iInt + 1;
  end
end

% Zero field interaction terms (S*D*S)
%--------------------------------------------------------------------------
if (nZFS>0)
  for iSpin = 1:nElSpins
    if System.S(iSpin)<1, continue; end
    % Get full 3x3 D matrix from spin system
    if System.fullD
      D_ = System.D(3*(iSpin-1)+(1:3),:);
    else
      D_ = diag(System.D(iSpin,:));
    end
    D_ = D_*1e6; % MHz -> Hz
    % Apply rotation if DFrame is given
    if any(System.DFrame(iSpin,:))
      R_M2D = erot(System.DFrame(iSpin,:)); % mol frame -> D frame
      R_D2M = R_M2D.'; % D frame -> mol frame
      D_ = R_D2M*D_*R_D2M.';
    end
    % Compute ISTO components
    if generalLiouvillian
      S_ = SpinOps(iSpin,:);
      [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(S_,S_);
    end
    [F0(iInt),F1(iInt,:),F2(iInt,:)] = tensor_cart2sph(D_);
    isFieldDep(iInt) = false;
    iInt = iInt + 1;
  end
end

% Electron-electron interaction terms (S1*J*S2)
%--------------------------------------------------------------------------
if nElSpins>1
  iCoupling = 1;
  for iEl1 = 1:nElSpins
    for iEl2 = iEl1+1:nElSpins
      % Construct matrix representing coupling tensor
      if System.fullee
        J_ = System.ee(3*(iCoupling-1)+(1:3),:); % MHz
      else
        J_ = diag(System.ee(iCoupling,:)); % MHz
      end
      if any(System.eeFrame(iCoupling,:))
        R_M2ee = erot(System.eeFrame(iCoupling,:)); % mol frame -> ee frame
        R_ee2M = R_M2ee.';  % ee frame -> mol frame
        J_ = R_ee2M*diag(System.ee(iCoupling,:))*R_ee2M.';
      end
      if generalLiouvillian
        [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(SpinOps(iEl1,:),SpinOps(iEl2,:));
      end
      [F0(iInt),F1(iInt,:),F2(iInt,:)] = tensor_cart2sph(J_*1e6); % MHz -> Hz
      isFieldDep(iInt) = false;
      iInt = iInt + 1;
      iCoupling = iCoupling + 1;
    end
  end
end

% Nuclear Zeeman interaction terms (-muN*B*gn*I/h)
%--------------------------------------------------------------------------
if generalLiouvillian
  if includeNuclearZeeman
    for iNucSpin = 1:nNucSpins
        I_ = SpinOps(nElSpins+iNucSpin,:);
        gn_ = System.gn(iNucSpin);
        [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(B0,I_);
        [F0(iInt),F1(iInt,:),F2(iInt,:)] = tensor_cart2sph(-gn_*nmagn/planck);
        isFieldDep(iInt) = true;
        iInt = iInt + 1;
    end
  end
else
  if includeNuclearZeeman
    gn0 = zeros(nNucSpins,1);
    for iNucSpin = 1:nNucSpins
      gn0(iNucSpin) = tensor_cart2sph(System.gn(iNucSpin));
      System.NZ0(iNucSpin) = nmagn*B0{3}*gn0(iNucSpin)/planck*2*pi; % -> angular freq.
    end
  else
    System.NZ0 = zeros(1,nNucSpins);
  end

  % Adaption for two nuclei, to feed to chili_liouvmatrix2
  if (System.nNuclei>=2)
    System.Ib = System.I(2);
    System.NZ0b = System.NZ0(2);
    System.HF0b = System.HF0(2);
    System.HF2b = System.HF2(:,2);
  end
end

% Nuclear electric quadrupole interaction terms
%--------------------------------------------------------------------------
if includeNuclearQuadrupole
  for iNucSpin = 1:nNucSpins
    Q_ = System.Q(iNucSpin,:)*1e6; % MHz -> Hz
    if any(System.QFrame(iNucSpin,:))
      R_M2Q = erot(System.QFrame(iNucSpin,:)); % mol frame -> Q frame
      R_Q2M = R_M2Q.'; % Q frame -> mol frame
      Q_ = R_Q2M*diag(Q_)*R_Q2M.';
    end
    if generalLiouvillian
      I_ = SpinOps(nElSpins+iNucSpin,:);
      [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(I_,I_);
    end
    [F0(iInt),F1(iInt,:),F2(iInt,:)] = tensor_cart2sph(Q_);
    isFieldDep(iInt) = false;
    iInt = iInt + 1;
  end
end

% Symmmetry tests
%--------------------------------------------------------------------------
% no need for including System.nNuclei==0 here since all interaction tensor
% components will be parsed, regardless of the number of nuclei
Symmetry.nobetatilts = allclose(F2(:,[2 4]),0);
Symmetry.tensorsCollinear = all(isreal(F2));
Symmetry.axialSystem = allclose(F2(:,1),0);

% Output
%--------------------------------------------------------------------------

% Pack up ISTOs
T.T0 = T0;
T.T1 = T1;
T.T2 = T2;

% Pack up ISTs
F.F0 = F0;
F.F1 = F1;
F.F2 = F2;

return
