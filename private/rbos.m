% rbos takes as input the Wigner D matrix, and the ISTO components and coefficients
% for each interaction term in the Hamiltonian. The output is 1 isotropic
% RBO and 25 anistropic rank-2 RBOs. The anisotropic RBOs are arranged as
% matrices in a 5-by-5 cell array in which the rows are the M components
% and the columns are the K components, in decreasing order from 2 to -2

function [Q0,Q1,Q2] = rbos(D1,D2,T0,T1,T2,F0,F1,F2)

nTerms = numel(F0);
dim = size(T0{1},1);

% construct the rank-0 isotropic RBO

Q0 = 0;
for iTerm = 1:nTerms
  Q0 = Q0 + conj(F0(iTerm))*T0{iTerm};
end


% construct the 9 rank-1 asymmetric RBOs

Q1 = cell(3,3);
for m_p = 1:3
  for m_q = 1:3
    
    Q1pq_ = 0;
    for m = 1:3
      for iTerm = 1:nTerms
        Q1pq_ = Q1pq_ + D1(m,m_p)*(-1)*conj(F1(iTerm,m_q))*T1{iTerm,m};
      end
    end
    Q1{m_p,m_q} = Q1pq_;
    
  end
end

        
% construct the 25 rank-2 anisotropic RBOs

Q2 = cell(5,5);
for m_p = 1:5
  for m_q = 1:5
    
    Q2pq_ = 0;
    for m = 1:5
      for iTerm = 1:nTerms
        Q2pq_ = Q2pq_ + D2(m,m_p)*conj(F2(iTerm,m_q))*T2{iTerm,m};
      end
    end
    Q2{m_p,m_q} = Q2pq_;
    
  end
end

% express the RBOs in Liouville space


LiouvilleSpace = true;
if LiouvilleSpace
  I = eye(dim);
  kronkron = @(A) kron(I,A)-kron(A.',I);
  %kronkron = @(A) kron(A,I)-kron(I,A.');
  Q0 = kronkron(Q0);
  for im1 = 1:3
    for im2 = 1:3
      Q1{im1,im2} = sparse(kronkron(Q1{im1,im2}));
    end
  end
  for im1 = 1:5
    for im2 = 1:5
      Q2{im1,im2} = sparse(kronkron(Q2{im1,im2}));
    end
  end
end
%{
idx = 1;
for i = 1:5
  for j = 1:5
    subplot(5,5,idx); 
    spy(Q2{i,j});
    xlabel(gca,sprintf('(%d,%d)',i,j));
    idx = idx + 1;
  end
end
%}
return