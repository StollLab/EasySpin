% rbos takes as input the ISTO components and coefficients for each
% interaction term in the Hamiltonian as well as the Euler angles describing
% the orientation of the spin system.
%
% The output is 1 isotropic
% RBO and 25 anistropic rank-2 RBOs. The anisotropic RBOs are arranged as
% matrices in a 5-by-5 cell array in which the rows are the M components
% and the columns are the K components, in decreasing order from 2 to -2

function [Q0G,Q1G,Q2G,Q0F,Q1F,Q2F] = rbos(T,F,angles,isFieldDep)

% Calculate Wigner rotation matrices
D1 = wignerd(1,angles(1),angles(2),angles(3));
D2 = wignerd(2,angles(1),angles(2),angles(3));

% Components of spin spherical tensor operators
T0 = T.T0;
T1 = T.T1;
T2 = T.T2;

% Components of interaction spherical tensors
F0 = F.F0;
F1 = F.F1;
F2 = F.F2;

nTerms = numel(F0);

% construct the rank-0 isotropic RBO
%---------------------------------------------------------------------
Q0G = sparse(0);
Q0F = sparse(0);
for iTerm = 1:nTerms
  if isFieldDep(iTerm)
    Q0G = Q0G + conj(F0(iTerm))*T0{iTerm};
  else
    Q0F = Q0F + conj(F0(iTerm))*T0{iTerm};
  end
end

% construct the 9 rank-1 asymmetric RBOs
%---------------------------------------------------------------------
if any(F1(:))
  Q1G = cell(3,3);
  Q1F = cell(3,3);
  for mp = 1:3
    for mq = 1:3
      
      QBmpmq = sparse(0);
      QGmpmq = sparse(0);
      for m = 1:3
        for iTerm = 1:nTerms
          if isFieldDep(iTerm)
            QBmpmq = QBmpmq + D1(m,mp)*(-1)*conj(F1(iTerm,mq))*T1{iTerm,m};
          else
            QGmpmq = QGmpmq + D1(m,mp)*(-1)*conj(F1(iTerm,mq))*T1{iTerm,m};
          end
        end
      end
      Q1G{mp,mq} = QBmpmq;
      Q1F{mp,mq} = QGmpmq;
      
    end
  end
else
  Q1G = {};
  Q1F = {};
end

% construct the 25 rank-2 anisotropic RBOs
%---------------------------------------------------------------------
Q2G = cell(5,5);
Q2F = cell(5,5);
for mp = 1:5
  for mq = 1:5
    
    QBmpmq = sparse(0);
    QGmpmq = sparse(0);
    for m = 1:5
      for iTerm = 1:nTerms
        if isFieldDep(iTerm)
          QBmpmq = QBmpmq + D2(m,mp)*conj(F2(iTerm,mq))*T2{iTerm,m};
        else
          QGmpmq = QGmpmq + D2(m,mp)*conj(F2(iTerm,mq))*T2{iTerm,m};
        end
      end
    end
    Q2G{mp,mq} = QBmpmq;
    Q2F{mp,mq} = QGmpmq;
    
  end
end

% Make sure Q*F are defined even in the absence of field-independent interactions
noFieldIndepTerm = all(isFieldDep);
if noFieldIndepTerm
  
  % rank-0
  nDim = length(T0{1});
  Q0F = sparse(nDim,nDim);
  
  % rank-1
  if any(F1(:))
    for i = 1:3
      for j = 1:3
        Q1F{i,j} = sparse(nDim,nDim);
      end
    end
  end
  
  % rank-2
  for i=1:5
    for j = 1:5
      Q2F{i,j} = sparse(nDim,nDim);
    end
  end
  
end
  
% Convert all operator matrices from Hilbert to Liouville space
Q0F = hil2liouv(Q0F);
Q0G = hil2liouv(Q0G);
Q1F = hil2liouv(Q1F);
Q1G = hil2liouv(Q1G);
Q2F = hil2liouv(Q2F);
Q2G = hil2liouv(Q2G);

return
