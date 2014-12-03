function [T0,T1,T2,F0,F1,F2] = magint(System,Experiment)

Sys = System;
Exp = Experiment;


% count the number of interaction terms in the spin Hamiltonian
% -------------------------------------------------------------------------


nElSpins = Sys.nElectrons;
Sys.I = nucspin(Sys.Nucs);
nNucSpins = Sys.nNuclei;

nZeeman = nElSpins;
nHyperfine = nElSpins*nNucSpins;
nZFS = 0;
if (numel(Sys.D) > 1 && any(Sys.D)) || (numel(Sys.fullD) > 1 && any(Sys.fullD)) 
  for iElectron = 1:nElSpins
    if (Sys.S(iElectron) > 1/2)
      nZFS = nZFS + 1;
    end
  end
end
    
nInteractions = nZeeman + nHyperfine + nZFS;

T0 = cell(nInteractions,1);
T1 = cell(nInteractions,3);
T2 = cell(nInteractions,5);
F0 = zeros(nInteractions,1);
F1 = zeros(nInteractions,3);
F2 = zeros(nInteractions,5);

iInteraction = 1;


% make input compatible with JerSpin istotensor and istocoef
% -------------------------------------------------------------------------

if (Sys.fullg == 0)
  g = diag(Sys.g);
end
B = {0 0 Exp.Field/1e3}; % mT -> T
Spins = [Sys.S Sys.I];
SpinOps = spinop(Spins);

A = cell(nNucSpins,1);
idx = 1;
for iNucSpin = 1:nNucSpins
  if (Sys.fullA ~= 0)
    A{iNucSpin} = Sys.A(idx:idx+2,1:3)*1e6; % MHz -> Hz
    idx = idx + 3;
  else
    A{iNucSpin} = diag(Sys.A(idx,1:3))*1e6; % MHz -> Hz
    idx = idx + 1;
  end
end

if (any(Sys.fullD)) || (any(Sys.D))
  if Sys.fullD == 0
    zfs = diag(Sys.D)*1e6; % MHz -> Hz
  else
    zfs = Sys.fullD*1e6;   % MHz -> Hz
  end
else
  zfs = 0;
end


%--------------------------------------------------------------------------
% electron Zeeman interaction terms
%--------------------------------------------------------------------------

for iSpin = 1:nElSpins
  S_ = SpinOps(iSpin,:);
  
  [T0_,T1_,T2_] = istotensor(B,S_);
  T0{iInteraction}   = T0_;
  T1(iInteraction,:) = T1_;
  T2(iInteraction,:) = T2_;
  [F0_,F1_,F2_] = istocoeff(g*bmagn/planck);
  F0(iInteraction)   = F0_;
  F1(iInteraction,:) = F1_;
  F2(iInteraction,:) = F2_;
  
  iInteraction = iInteraction + 1;
end

%--------------------------------------------------------------------------
% hyperfine interaction terms
%--------------------------------------------------------------------------

nNuc = 1;
for iSpin1 = 1:nElSpins
  S_ = SpinOps(iSpin1,:);
  for iSpin2 = nElSpins+1:numel(Spins)
    I_ = SpinOps(iSpin2,:);
    
    [T0_,T1_,T2_] = istotensor(S_,I_);
    T0{iInteraction}   = T0_;
    T1(iInteraction,:) = T1_;
    T2(iInteraction,:) = T2_;
    
    [F0_,F1_,F2_] = istocoeff(A{nNuc});
    F0(iInteraction) =   F0_;
    F1(iInteraction,:) = F1_;
    F2(iInteraction,:) = F2_;
    
    nNuc = nNuc + 1;
    iInteraction = iInteraction + 1;
  end
  nNuc = 1;
end

%--------------------------------------------------------------------------
% zero field splitting
%--------------------------------------------------------------------------

if nZFS > 0
  for iSpin = 1:nElSpins
    if Spins(iSpin) > 1/2
      S_ = SpinOps(iSpin,:);
      
      [T0_,T1_,T2_] = istotensor(S_,S_);
      T0{iInteraction}   = T0_;
      T1(iInteraction,:) = T1_;
      T2(iInteraction,:) = T2_;
      [F0_,F1_,F2_] = istocoeff(zfs);
      F0(iInteraction)   = F0_;
      F1(iInteraction,:) = F1_;
      F2(iInteraction,:) = F2_;
      
      iInteraction = iInteraction + 1;
    end
  end
end


return