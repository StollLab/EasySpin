% Calculate Hamiltonian derivatives with respect to D and E, taking into
% account correlation if present; premultipy with linewidth.

% Required input fields
%   Sys.DStrain:
%     nElectrons x 2 array with [FWHM_D FWHM_E] on each row.
%   Sys.DStrainCorr:
%     Array with nElectron correlation coefficient between D and E (between -1 and +1).
%   Sys.DFrame:
%     Euler angles describing D tensor orientation in molecular frame

function dHdDE = strains_de(Sys)

% Return if no D/E strain present
if ~any(Sys.DStrain(:))
  dHdDE = {};
  return
end

% Loop over all electron spins
nElectrons = size(Sys.DStrain,1);
dHdDE = cell(1,2*nElectrons);  % preallocate maximal possible size
idx = 0;
for iEl = 1:nElectrons
  
  % Skip if no D strain is given for this electron spin
  if ~any(Sys.DStrain(iEl,1:2))
    continue
  end
  
  % Strain information
  Delta = Sys.DStrain(iEl,:);  % FWHMs for D and E
  rDE = Sys.DStrainCorr(iEl);  % correlation coefficient between D and E
  if rDE<-1 || rDE>1
    error('Electron spin %d: D/E correlation coefficient must be between -1 and 1, verify your input',iEl);
  end

  % Construct spin operators in molecular frame
  SxM = sop(Sys,[iEl,1]);
  SyM = sop(Sys,[iEl,2]);
  SzM = sop(Sys,[iEl,3]);
  
  % Transform spin operators to D frame
  if any(Sys.DFrame(iEl,:))
    R_M2D = erot(Sys.DFrame(iEl,:));  % molecular frame -> D frame  -----------------------M2D pr D2M? confused due to g2M in strains_ga
    SxD = R_M2D(1,1)*SxM + R_M2D(1,2)*SyM + R_M2D(1,3)*SzM;
    SyD = R_M2D(2,1)*SxM + R_M2D(2,2)*SyM + R_M2D(2,3)*SzM;
    SzD = R_M2D(3,1)*SxM + R_M2D(3,2)*SyM + R_M2D(3,3)*SzM;
  else
    % D frame aligns with molecular frame
    SxD = SxM;
    SyD = SyM;
    SzD = SzM;
  end
  
  % Calculate derivatives of Hamiltonian w.r.t. D and E
  % Spin Hamiltonian terms (in D tensor frame xD, yD, zD)
  %   H = D*(SzD^2-S(S+1)/3) + E*(SxD^2-SyD^2)
  %     = D*(2*SzD^2-SxD^2-SyD^2)/3 + E*(SxD^2-SyD^2)
  dHdD = (2*SzD^2 - SxD^2 - SyD^2)/3;
  dHdE = SxD^2 - SyD^2;
  
  % Diagonalize covariance matrix in the case of non-zero correlation
  if rDE~=0
    % Transform correlated D-E strain to uncorrelated coordinates
    logmsg(1,'  correlated D strain for electron spin %d (r = %f)',iEl,rDE);
    % Construct covariance matrix
    R12 = rDE*Delta(1)*Delta(2);  % covariance
    covMatrix = [Delta(1)^2 R12; R12 Delta(2)^2];
    % Diagonalize covariance matrix, get derivatives along principal directions
    [V,sigma2] = eig(covMatrix);
    Delta = sqrt(diag(sigma2));
    dHdDE1 = Delta(1)*(V(1,1)*dHdD + V(1,2)*dHdE);
    dHdDE2 = Delta(2)*(V(2,1)*dHdD + V(2,2)*dHdE);
  else
    dHdDE1 = Delta(1)*dHdD;
    dHdDE2 = Delta(2)*dHdE;
  end

  % Store in list
  idx = idx+1;
  dHdDE{idx} = dHdDE1;
  idx = idx+1;
  dHdDE{idx} = dHdDE2;
  
end

dHdDE = dHdDE(1:idx);

end
