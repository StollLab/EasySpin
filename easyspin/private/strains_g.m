% Calculate Hamiltonian derivatives with respect to the g Principle Axis values, taking into
% account correlation if present; premultipy with linewidth.

% Required input fields
%   Sys.gStrain:
%     nElectrons x 3 array with [gxx gyy gzz] on each row.
%   Sys.gStrainCorr:
%     Array with nElectron correlation coefficient between gxx, gyy and gzz (between -1 and +1).
%   Sys.gFrame:
%     Euler angles describing g tensor orientation in molecular frame

function dHdDg = strains_g(Sys)

% Return if no g strain present
if ~any(Sys.gStrain(:))
  dHdDg = {};
  return
end

% Loop over all electron spins
nElectrons = size(Sys.gStrain,1);
dHdDg = cell(3,nElectrons);  % preallocate maximal possible size 3*number of electrons
idx = 0;
for iEl = 1:nElectrons
  
  % Skip if no g strain is given for this electron spin
  if ~any(Sys.gStrain(iEl,:))
    continue
  end
  
  % Strain information
  Delta = Sys.gStrain(iEl,:);  % FWHMs for g
  rgg = Sys.gStrainCorr(:,:,iEl);  % correlation coefficient between g components (3x3 matrix per electron) --> should be a 3D matrix?
  if find(rgg)<-1 || find(rgg)>1
    error('Electron spin %d: gi/gj correlation coefficient must be between -1 and 1, verify your input',iEl);
  end

  % Construct spin operators in molecular frame
  SxM = sop(Sys,[iEl,1]);
  SyM = sop(Sys,[iEl,2]);
  SzM = sop(Sys,[iEl,3]);
  
  % Transform spin operators to g frame
  if any(Sys.gFrame(iEl,:))
    R_g2M = erot(Sys.gFrame(iEl,:)).';  % molecular frame -> g Frame


    SxD = R_M2D(1,1)*SxM + R_M2D(1,2)*SyM + R_M2D(1,3)*SzM;
    SyD = R_M2D(2,1)*SxM + R_M2D(2,2)*SyM + R_M2D(2,3)*SzM;
    SzD = R_M2D(3,1)*SxM + R_M2D(3,2)*SyM + R_M2D(3,3)*SzM;
  else
    % g frame aligns with molecular frame
    SxD = SxM;
    SyD = SyM;
    SzD = SzM;
  end
  
  % Calculate derivatives of Hamiltonian w.r.t. gxx, gyy and gzz
  
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
% Calculates quantities needed by resfields for g/A strain broadening

function [usegStrain,useAStrain,simplegStrain,lw2_gA,ops] = strains_ga(CoreSys,Exp,mwFreq,Transitions,nTransitions,kH0,kmuzM)

usegStrain = any(CoreSys.gStrain(:));
useAStrain = CoreSys.nNuclei>0 && any(CoreSys.AStrain);

ops = struct;
lw2_gA = [];

% g-A strain
%-------------------------------------------------
% g strain tensor is taken to be aligned with the g tensor
% A strain tensor is taken to be aligned with the A tensor
% g strain can be specified for each electron spin
% A strain is limited to the first electron and first nuclear spin
simplegStrain = CoreSys.nElectrons==1;
if usegStrain
  logmsg(1,'  g strain present');
  lw_g = mwFreq*CoreSys.gStrain./CoreSys.g;  % MHz
  for Ne = CoreSys.nElectrons:-1:1
    if any(CoreSys.gFrame(Ne,:))
      R_g2M{Ne} = erot(CoreSys.gFrame(Ne,:)).';  % g frame -> molecular frame [question about why .' not consistent with strains_de]
    else
      R_g2M{Ne} = eye(3);
    end
  end
  if ~simplegStrain
    logmsg(1,'  multiple g strains present');
    for Ne = CoreSys.nElectrons:-1:1
      kSxM{Ne} = sop(CoreSys,[Ne,1]);
      kSyM{Ne} = sop(CoreSys,[Ne,2]);
      kSzM{Ne} = sop(CoreSys,[Ne,3]);
    end
    ops.kSxM = kSxM;
    ops.kSyM = kSyM;
    ops.kSzM = kSzM;
  end
else
  lw_g = zeros(CoreSys.nElectrons,3);
end

if useAStrain
  if usegStrain
    if any(CoreSys.gFrame(1,:)~=CoreSys.AFrame(1,1:3))
      error('For g/A strain, Sys.gFrame and Sys.AFrame must be collinear.');
    end
    R_strain2M = R_g2M;
  else
    R_A2M = erot(CoreSys.AFrame).';  % strain tensor frame -> mol. frame
    R_strain2M = R_A2M;
  end
  lw_A = diag(CoreSys.AStrain);
  % Diagonalize Hamiltonian at center field.
  centerB = mean(Exp.Range);
  [Vecs,E] = eig(kH0 - centerB*kmuzM);
  [~,idx] = sort(real(diag(E)));
  Vecs = Vecs(:,idx);
  % Calculate effective mI of nucleus 1 for all eigenstates.
  mI1 = real(diag(Vecs'*sop(CoreSys,[2,3])*Vecs));
  mI1Tr = mean(mI1(Transitions),2);
  % compute A strain array
  lw_A = reshape(mI1Tr(:,ones(1,9)).',[3,3,nTransitions]).*...
    repmat(lw_A,[1,1,nTransitions]);

  % Combine g and A strain
  rho = CoreSys.gAStrainCorr;  % correlation coefficient
  for Ne = CoreSys.nElectrons:-1:1
    lw2_gA{Ne} = repmat(diag(lw_g(Ne,:).^2),[1,1,nTransitions]) + lw_A.^2 + ...
      2*rho*repmat(diag(lw_g(Ne,:)),[1,1,nTransitions]).*lw_A;
    for tr = 1:nTransitions
      lw2_gA{Ne}(:,:,tr) = R_strain2M{Ne}*lw2_gA{Ne}(:,:,tr)*R_strain2M{Ne}.';
    end
  end
else
  if usegStrain
    for Ne = CoreSys.nElectrons:-1:1
      lw2_gA{Ne} = repmat(R_g2M{Ne}*diag(lw_g(Ne,:).^2)*R_g2M{Ne}.',[1,1,nTransitions]);
    end
  end
end
% lw2_gA = a (cell array of) 3D array with 3x3 strain line-width matrices
% for each transition stacked along the third dimension.
