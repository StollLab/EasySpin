% rbos takes as input the Wigner D matrix, and the ISTO components and coefficients
% for each interaction term in the Hamiltonian. The output is 1 isotropic
% RBO and 25 anistropic rank-2 RBOs. The anisotropic RBOs are arranged as
% matrices in a 5-by-5 cell array in which the rows are the M components
% and the columns are the K components, in decreasing order from 2 to -2

function [Q0B,Q1B,Q2B,Q0G,Q1G,Q2G] = rbos(D1,D2,T0,T1,T2,F0,F1,F2,isFieldDep)

nTerms = numel(F0);

% construct the rank-0 isotropic RBO
%---------------------------------------------------------------------
Q0B = sparse(0);
Q0G = sparse(0);
for iTerm = 1:nTerms
  if isFieldDep(iTerm)
    Q0B = Q0B + conj(F0(iTerm))*T0{iTerm};
  else
    Q0G = Q0G + conj(F0(iTerm))*T0{iTerm};
  end
end

% construct the 9 rank-1 asymmetric RBOs
%---------------------------------------------------------------------
if any(F1(:))
  Q1B = cell(3,3);
  Q1G = cell(3,3);
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
      Q1B{mp,mq} = QBmpmq;
      Q1G{mp,mq} = QGmpmq;
      
    end
  end
else
  Q1B = {};
  Q1G = {};
end

% construct the 25 rank-2 anisotropic RBOs
%---------------------------------------------------------------------
Q2B = cell(5,5);
Q2G = cell(5,5);
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
    Q2B{mp,mq} = QBmpmq;
    Q2G{mp,mq} = QGmpmq;
    
  end
end

% check for existence of any field-independent interactions
if ~any(isFieldDep==0)
  
  % rank-0
  nDim = length(T0{1});
  Q0G = sparse(nDim,nDim);
  
  % rank-1
  if any(F1(:))
    for i = 1:3
      for j = 1:3
        Q1G{i,j} = sparse(nDim,nDim);
      end
    end
  end
  
  % rank-2
  for i=1:5
    for j = 1:5
      Q2G{i,j} = sparse(nDim,nDim);
    end
  end
  
end
  


return
