% strains_ops  Strain mode operators
%
%   [G0,Gmu] = strains_ops(Sys,CoreSys,coreNuclei,sparseFlag)
%
%   Constructs, for each strain mode k, the derivative of the spin Hamiltonian
%   along the mode, G_k = sum_i Q(i,k)*dH/dp_i = G0{k} - B*zL.'*[Gmu{:,k}].
%
%   Input:
%     Sys          validated spin system with Sys.StrainData
%     CoreSys      spin system whose spin space the operators are built on
%     coreNuclei   indices (in Sys) of the nuclei that are kept in CoreSys
%     sparseFlag   true for sparse output matrices
%
%   Output:
%     G0           1 x m cell array of field-independent operators (MHz)
%     Gmu          3 x m cell array of magnetic-moment operators (MHz/mT),
%                  molecular-frame components x, y, z

function [G0,Gmu] = strains_ops(Sys,CoreSys,coreNuclei,sparseFlag)

Tensors = Sys.StrainData.Tensors;
nModes = size(Sys.StrainData.Q,2);
nStates = hsdim(CoreSys);
nEl = CoreSys.nElectrons;

% Map nucleus indices from Sys to CoreSys
nucMap = zeros(1,Sys.nNuclei);
nucMap(coreNuclei) = 1:numel(coreNuclei);

G0 = cell(1,nModes);
Gmu = cell(3,nModes);
for k = nModes:-1:1
  G0{k} = sparse(nStates,nStates);
  for a = 3:-1:1
    Gmu{a,k} = sparse(nStates,nStates);
  end
end

pre_e = -bmagn/planck/1e9; % MHz/mT
for iTensor = 1:numel(Tensors)
  d = Tensors(iTensor);
  idx = d.idx;

  % Spin indices in CoreSys
  switch d.type
    case {'g','D'}
      spins = idx;
    case {'Q','sigma'}
      spins = nEl + corenucleus(idx);
    case 'A'
      spins = [idx(1), nEl+corenucleus(idx(2))];
    case 'ee'
      spins = idx;
    case 'nn'
      spins = nEl + [corenucleus(idx(1)) corenucleus(idx(2))];
  end

  % Spin operator components
  for c = 3:-1:1
    S1{c} = sop(CoreSys,[spins(1),c],'sparse');
    S2{c} = sop(CoreSys,[spins(end),c],'sparse');
  end

  for k = 1:nModes
    % Derivative of the tensor along mode k
    M = d.M(:,:,k);
    if ~any(M(:)), continue; end

    switch d.type
      case 'g'
        for a = 1:3
          for j = 1:3
            Gmu{a,k} = Gmu{a,k} + pre_e*M(a,j)*S1{j};
          end
        end
      case 'sigma'
        n = idx;
        pre_n = nmagn/planck/1e9*Sys.gn(n)*Sys.gnscale(n); % MHz/mT
        for a = 1:3
          for j = 1:3
            Gmu{a,k} = Gmu{a,k} + pre_n*M(a,j)*S1{j};
          end
        end
      otherwise % A, D, Q, ee, nn: bilinear terms S1*M*S2
        for i = 1:3
          for j = 1:3
            G0{k} = G0{k} + M(i,j)*(S1{i}*S2{j});
          end
        end
    end
  end
end

if ~sparseFlag
  G0 = cellfun(@full,G0,'UniformOutput',false);
  Gmu = cellfun(@full,Gmu,'UniformOutput',false);
end

  %-----------------------------------------------------------------------------
  function nc = corenucleus(n)
  nc = nucMap(n);
  if nc==0
    error('Strain parameter refers to nucleus %d, which is not in the core system.',n);
  end
  end

end
