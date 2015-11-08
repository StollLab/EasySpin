% rbos takes as input the Wigner D matrix, and the ISTO components and coefficients
% for each interaction term in the Hamiltonian. The output is 1 isotropic
% RBO and 25 anistropic rank-2 RBOs. The anisotropic RBOs are arranged as
% matrices in a 5-by-5 cell array in which the rows are the M components
% and the columns are the K components, in decreasing order from 2 to -2

function [Q0,Q1,Q2] = rbos(D1,D2,T0,T1,T2,F0,F1,F2)

nTerms = numel(F0);

% construct the rank-0 isotropic RBO
%---------------------------------------------------------------------
Q0 = sparse(0);
for iTerm = 1:nTerms
  Q0 = Q0 + conj(F0(iTerm))*T0{iTerm};
end

% construct the 9 rank-1 asymmetric RBOs
%---------------------------------------------------------------------
if any(F1(:))
  Q1 = cell(3,3);
  for mp = 1:3
    for mq = 1:3
      
      Qmpmq = sparse(0);
      for m = 1:3
        for iTerm = 1:nTerms
          Qmpmq = Qmpmq + D1(m,mp)*(-1)*conj(F1(iTerm,mq))*T1{iTerm,m};
        end
      end
      Q1{mp,mq} = Qmpmq;
      
    end
  end
else
  Q1 = {};
end

% construct the 25 rank-2 anisotropic RBOs
%---------------------------------------------------------------------
Q2 = cell(5,5);
for mp = 1:5
  for mq = 1:5
    
    Qmpmq = sparse(0);
    for m = 1:5
      for iTerm = 1:nTerms
        Qmpmq = Qmpmq + D2(m,mp)*conj(F2(iTerm,mq))*T2{iTerm,m};
      end
    end
    Q2{mp,mq} = Qmpmq;
    
  end
end

% express the RBOs in Liouville space
%---------------------------------------------------------------------
nI = length(Q0);
kronkron = @(A) speyekron(nI,A)-spkroneye(A.',nI);
Q0 = kronkron(Q0);
if ~isempty(Q1)
  for im1 = 1:3
    for im2 = 1:3
      Q1{im1,im2} = kronkron(Q1{im1,im2});
    end
  end
end
for im1 = 1:5
  for im2 = 1:5
    Q2{im1,im2} = kronkron(Q2{im1,im2});
  end
end

return
