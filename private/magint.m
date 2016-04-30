function [T0,T1,T2,F0,F1,F2] = magint(System,SpinOps,CenterField,IncludeNuclearZeeman)

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

if IncludeNuclearZeeman
  nInteractions = nInteractions + nNucSpins;
end

IncludeNuclearQuadrupole = false;
if IncludeNuclearQuadrupole
  nInteractions = nInteractions + nNucSpins;
end

T0 = cell(nInteractions,1);
T1 = cell(nInteractions,3);
T2 = cell(nInteractions,5);
F0 = zeros(nInteractions,1);
F1 = zeros(nInteractions,3);
F2 = zeros(nInteractions,5);

iInt = 1;

B0 = {0 0 CenterField/1e3}; % mT -> T

% Electron Zeeman interaction terms (muB*B*g*S/h)
%--------------------------------------------------------------------------
for iElSpin = 1:nElSpins
  if ~System.fullg
    g = diag(System.g(iElSpin,:));
  else
    g = System.g((iElSpin-1)*3+(1:3),:);
  end
  [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(B0,SpinOps(iElSpin,:));
  [F0(iInt),F1(iInt,:),F2(iInt,:)] = istocoeff(g*bmagn/planck); % Hz
  iInt = iInt + 1;
end

% Hyperfine interaction terms (S*A*I)
%--------------------------------------------------------------------------
for iElSpin = 1:nElSpins
  eidx = 3*(iElSpin-1)+(1:3);
  S_ = SpinOps(iElSpin,:);
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
    I_ = SpinOps(nElSpins+iNucSpin,:);
    [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(S_,I_);
    [F0(iInt),F1(iInt,:),F2(iInt,:)] = istocoeff(A_);
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
    S_ = SpinOps(iSpin,:);
    [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(S_,S_);
    [F0(iInt),F1(iInt,:),F2(iInt,:)] = istocoeff(D_);
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
      [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(SpinOps(iEl1,:),SpinOps(iEl2,:));
      [F0(iInt),F1(iInt,:),F2(iInt,:)] = istocoeff(J_*1e6); % MHz -> Hz
      iInt = iInt + 1;
      iCoupling = iCoupling + 1;
    end
  end
end

% Nuclear Zeeman interaction terms (-muN*B*gn*I/h)
%--------------------------------------------------------------------------
if IncludeNuclearZeeman
  for iNucSpin = 1:nNucSpins
    I_ = SpinOps(nElSpins+iNucSpin,:);
    gn_ = System.gn(iNucSpin);
    [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(B0,I_);
    [F0(iInt),F1(iInt,:),F2(iInt,:)] = istocoeff(-gn_*nmagn/planck);
    iInt = iInt + 1;
  end
end

% Nuclear electric quadrupole interaction terms
%--------------------------------------------------------------------------
if IncludeNuclearQuadrupole
  for iNucSpin = 1:nNucSpins
    I_ = SpinOps(nElSpins+iNucSpin,:);
    Q_ = System.Q(iNucSpin,:)*1e6; % MHz -> Hz
    if any(System.QFrame(iNucSpin,:))
      R_M2Q = erot(System.QFrame(iNucSpin,:)); % mol frame -> Q frame
      R_Q2M = R_M2Q.'; % Q frame -> mol frame
      Q_ = R_Q2M*diag(Q_)*R_Q2M.';
    end
    [T0{iInt},T1(iInt,:),T2(iInt,:)] = istotensor(I_,I_);
    [F0(iInt),F1(iInt,:),F2(iInt,:)] = istocoeff(Q_);
    iInt = iInt + 1;
  end
end

return
