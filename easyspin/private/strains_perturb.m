% strains_perturb  Strain widths for perturbation-theory line positions
%
%   dE2 = strains_perturb(Sys,nB,E0,EZ,secondOrder)
%
%   Computes, for each allowed transition and orientation, the sum over strain
%   modes of the squared derivatives of the transition energy at fixed field,
%   using analytical derivatives of the first- or second-order perturbation
%   expressions of Iwasaki, J. Magn. Reson. 16, 417 (1974), Eqs. [20]-[28],
%   https://doi.org/10.1016/0022-2364(74)90223-6
%
%   Input:
%     Sys          validated spin system with one electron spin and
%                  Sys.StrainData
%     nB           3 x nOri array of field directions, molecular frame
%     E0           electron Zeeman energy used in the second-order terms, MHz;
%                  scalar (field sweeps) or 1 x nOri (frequency sweeps)
%     EZ           nRows x nOri array, coefficient of dgeff/geff in the energy
%                  derivative, MHz
%                    field sweeps:      E0 - (all corrections)
%                    frequency sweeps:  E0 - (second-order corrections)
%     secondOrder  true for second-order perturbation theory
%
%   Output:
%     dE2          nRows x nOri array of sum_k dE_k^2, MHz^2
%
%   Rows are ordered by mS <-> mS-1 transition (mS = S, S-1, ..., -S+1), and
%   within each by nuclear sublevel combination as given by allcombinations.

function dE2 = strains_perturb(Sys,nB,E0,EZ,secondOrder)

S = Sys.S;
highSpin = S>1/2;
nNuclei = Sys.nNuclei;
nOri = size(nB,2);
if isscalar(E0), E0 = E0*ones(1,nOri); end

% Tensors in the molecular frame, as in resfields_perturb and resfreqs_perturb
if Sys.fullg
  g = Sys.g;
else
  R_g2M = erot(Sys.gFrame).'; % g frame -> molecular frame
  g = R_g2M*diag(Sys.g)*R_g2M.';
end
if highSpin
  if Sys.fullD
    D = Sys.D;
  else
    R_D2M = erot(Sys.DFrame).'; % D frame -> molecular frame
    D = R_D2M*diag(Sys.D)*R_D2M.';
  end
  D = D - eye(3)*trace(D)/3;
end
I = Sys.I;
for n = nNuclei:-1:1
  if Sys.fullA
    % Iwasaki's Hamiltonian is I.A.S, EasySpin's is S.A.I, so transpose
    A{n} = Sys.A((n-1)*3+(1:3),:).';
  else
    R_A2M = erot(Sys.AFrame(n,:)).'; % A frame -> molecular frame
    A{n} = R_A2M*diag(Sys.A(n,:))*R_A2M.';
  end
  detA(n) = det(A{n});
  invA{n} = inv(A{n});
  mI{n} = -I(n):I(n);
end
if nNuclei>0
  mIc = allcombinations(mI{:});
  II1 = I.*(I+1);
else
  mIc = zeros(1,0);
end
nNucSublevels = size(mIc,1);

% Tensor derivatives along the strain modes
nModes = size(Sys.StrainData.Q,2);
Mg = zeros(3,3,nModes);
MD = zeros(3,3,nModes);
MA = repmat({zeros(3,3,nModes)},1,nNuclei);
for t = Sys.StrainData.Tensors
  switch t.type
    case 'g'
      Mg = t.M;
    case 'D'
      MD = t.M;
      for k = 1:nModes
        MD(:,:,k) = MD(:,:,k) - eye(3)*trace(MD(:,:,k))/3; % D is made traceless
      end
    case 'A'
      MA{t.idx(2)} = permute(t.M,[2 1 3]); % transpose, as for A
    otherwise
      error('Strains of %s are not supported by perturbation theory.',t.type);
  end
end
gStrained = any(Mg(:));
AStrained = cellfun(@(M)any(M(:)),MA);

dE2 = zeros(2*S*nNucSublevels,nOri);
for iOri = 1:nOri
  n0 = nB(:,iOri);
  v = g.'*n0;
  geff = norm(v);
  u = v/geff;

  % Derivatives of geff and of the quantization axis u
  dgeff = zeros(1,nModes);
  du = zeros(3,nModes);
  if gStrained
    for k = 1:nModes
      dv = Mg(:,:,k).'*n0;
      dgeff(k) = u.'*dv;
      du(:,k) = (dv - u*dgeff(k))/geff;
    end
  end

  % Zero-field splitting
  if highSpin
    Du = D*u;
    uDu = u.'*Du;
    dDu = zeros(3,nModes);
    duDu = zeros(1,nModes);
    dtrDD = zeros(1,nModes);
    for k = 1:nModes
      dDu(:,k) = MD(:,:,k)*u + D*du(:,k);
      duDu(k) = u.'*MD(:,:,k)*u + u.'*(D+D.')*du(:,k);
      dtrDD(k) = 2*trace(D*MD(:,:,k));
    end
    duDDu = 2*Du.'*dDu;
    dD1sq = duDDu - 2*uDu*duDu;
    dD2sq = 2*dtrDD + 2*uDu*duDu - 4*duDDu;
  end

  % Hyperfine couplings, one row per nucleus
  dnK = zeros(nNuclei,nModes);
  dA1sq = zeros(nNuclei,nModes);
  dA2 = zeros(nNuclei,nModes);
  dA3 = zeros(nNuclei,nModes);
  dDA = zeros(nNuclei,nModes);
  for n = 1:nNuclei
    A_ = A{n};
    K = A_*u;
    nK = norm(K);
    kv = K/nK;
    Ak = A_.'*kv;
    if gStrained || AStrained(n)
      M = MA{n};
      dk = zeros(3,nModes);
      dAk = zeros(3,nModes);
      for k = 1:nModes
        a = M(:,:,k)*u + A_*du(:,k);
        dnK(n,k) = kv.'*a;
        dk(:,k) = (a - kv*dnK(n,k))/nK;
        dAk(:,k) = M(:,:,k).'*kv + A_.'*dk(:,k);
        dinvA = -invA{n}*M(:,:,k)*invA{n};
        ddetA = detA(n)*trace(invA{n}*M(:,:,k));
        dA2(n,k) = ddetA*(u.'*invA{n}*kv) + ...
          detA(n)*(du(:,k).'*invA{n}*kv + u.'*dinvA*kv + u.'*invA{n}*dk(:,k));
        dA3(n,k) = 2*trace(A_.'*M(:,:,k));
      end
      dkAAk = 2*Ak.'*dAk;
      dA1sq(n,:) = dkAAk - 2*nK*dnK(n,:); % kAu = nK
      dA3(n,:) = dA3(n,:) - dkAAk;
    else
      dAk = zeros(3,nModes);
    end
    if highSpin
      dDA(n,:) = Ak.'*dDu + Du.'*dAk - nK*duDu - uDu*dnK(n,:);
    end
  end

  % Loop over all mS <-> mS-1 transitions
  for imS = 1:2*S
    mS = S + 1 - imS;
    dE = zeros(nNucSublevels,nModes);

    % first order
    if highSpin
      dE = dE - duDu/2*(3-6*mS);
    end
    if nNuclei>0
      dE = dE + mIc*dnK;
    end

    % second order
    if secondOrder
      if highSpin
        c1 = 4*S*(S+1) - 3*(8*mS^2-8*mS+3);
        c2 = 2*S*(S+1) - 3*(2*mS^2-2*mS+1);
        dE = dE - (c1*dD1sq - c2/4*dD2sq)/(2*E0(iOri));
      end
      if nNuclei>0
        x = mIc.^2*dA1sq - (1-2*mS)*mIc*dA2 + (II1-mIc.^2)*dA3/2;
        dE = dE + x/(2*E0(iOri));
        if highSpin
          dE = dE - (3-6*mS)*mIc*dDA/E0(iOri);
        end
      end
    end

    % electron Zeeman
    rows = (imS-1)*nNucSublevels + (1:nNucSublevels);
    dE = dE + EZ(rows,iOri).*dgeff/geff;

    dE2(rows,iOri) = sum(dE.^2,2);
  end

end

end
