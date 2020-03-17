% Calculate Hamiltonian derivatives with respect to D and E, taking into
% account correlation if present), premultipy with strains.

% Sys.DStrain: nElectrons x 2 or nElectrons x 3 array
% with [FWHM_D FWHM_E rDE] on each row. rDE is the correlation coefficient
% between D and E.

function [useDStrain,dHdD,dHdE] = getdstrainops(Sys)

useDStrain = any(Sys.DStrain(:));
dHdD = cell(1,Sys.nElectrons);
dHdE = cell(1,Sys.nElectrons);

if ~useDStrain
  return
end

% Loop over all electron spins
for iEl = 1:Sys.nElectrons
  
  % Skip if no D strain is given for this electron spin
  if ~any(Sys.DStrain(iEl,1:2))
    dHdD{iEl} = 0;
    dHdE{iEl} = 0;
    continue
  end
  
  % Construct Zeeman operators in molecular frame
  SxM_ = sop(Sys,[iEl,1]);
  SyM_ = sop(Sys,[iEl,2]);
  SzM_ = sop(Sys,[iEl,3]);
  
  % Construct Zeeman operators in D frame, and square
  if any(Sys.DFrame(iEl,:))
    R_M2D = erot(Sys.DFrame(iEl,:)); % molecular frame -> D frame
    SxD2_ = (R_M2D(1,1)*SxM_ + R_M2D(1,2)*SyM_ + R_M2D(1,3)*SzM_)^2;
    SyD2_ = (R_M2D(2,1)*SxM_ + R_M2D(2,2)*SyM_ + R_M2D(2,3)*SzM_)^2;
    SzD2_ = (R_M2D(3,1)*SxM_ + R_M2D(3,2)*SyM_ + R_M2D(3,3)*SzM_)^2;
  else
    % D frame aligns with molecular frame
    SxD2_ = SxM_^2;
    SyD2_ = SyM_^2;
    SzD2_ = SzM_^2;
  end
  
  % Calculate derivatives of Hamiltonian w.r.t. D and E
  %   Spin Hamiltonian terms (in D tensor frame xD, yD, zD)
  %   H = D*(SzD^2-S(S+1)/3) + E*(SxD^2-SyD^2)
  %     = D*(2*SzD^2-SxD^2-SyD^2)/3 + E*(SxD^2-SyD^2)
  dHdD_ = (2*SzD2_-SxD2_-SyD2_)/3;
  dHdE_ = SxD2_-SyD2_;
  
  % Compute Hamiltonian derivatives, pre-multiply with strain FWHMs.
  DeltaD = Sys.DStrain(iEl,1);
  DeltaE = Sys.DStrain(iEl,2);
  rDE = Sys.DStrainCorr(iEl); % correlation coefficient between D and E
  if (rDE~=0)
    % Transform correlated D-E strain to uncorrelated coordinates
    logmsg(1,'  correlated D strain for electron spin %d (r = %f)',iEl,rDE);
    % Construct and diagonalize covariance matrix
    R12 = rDE*DeltaD*DeltaE;
    CovMatrix = [DeltaD^2 R12; R12 DeltaE^2];
    [V,L] = eig(CovMatrix);
    L = sqrt(diag(L));
    dHdD{iEl} = L(1)*(V(1,1)*dHdD_ + V(1,2)*dHdE_);
    dHdE{iEl} = L(2)*(V(2,1)*dHdD_ + V(2,2)*dHdE_);
  else
    dHdD{iEl} = DeltaD*dHdD_;
    dHdE{iEl} = DeltaE*dHdE_;
  end
  
end
