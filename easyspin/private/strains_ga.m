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
  for e = CoreSys.nElectrons:-1:1
    if any(CoreSys.gFrame(e,:))
      R_g2M{e} = erot(CoreSys.gFrame(e,:)).';  % g frame -> molecular frame
    else
      R_g2M{e} = eye(3);
    end
  end
  if ~simplegStrain
    logmsg(1,'  multiple g strains present');
    for e = CoreSys.nElectrons:-1:1
      kSxM{e} = sop(CoreSys,[e,1]);
      kSyM{e} = sop(CoreSys,[e,2]);
      kSzM{e} = sop(CoreSys,[e,3]);
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
  for e = CoreSys.nElectrons:-1:1
    lw2_gA{e} = repmat(diag(lw_g(e,:).^2),[1,1,nTransitions]) + lw_A.^2 + ...
      2*rho*repmat(diag(lw_g(e,:)),[1,1,nTransitions]).*lw_A;
    for tr = 1:nTransitions
      lw2_gA{e}(:,:,tr) = R_strain2M{e}*lw2_gA{e}(:,:,tr)*R_strain2M{e}.';
    end
  end
else
  if usegStrain
    for e = CoreSys.nElectrons:-1:1
      lw2_gA{e} = repmat(R_g2M{e}*diag(lw_g(e,:).^2)*R_g2M{e}.',[1,1,nTransitions]);
    end
  end
end
% lw2_gA = a (cell array of) 3D array with 3x3 strain line-width matrices
% for each transition stacked along the third dimension.
