% isot2stev  Transformation from ISTOs to Stevens operators
%   
%   C = isto2stev(k)
%
% Provide transformation matrix to transform a rank-k irreducible spherical tensor
% operator (ISTO) T_k to a vector of rank_k extended Steven operators O_k with
%
%   O_k = C *T_k
%
% T_k is the ISTO of rank k with 2*k+1 components T(k,q=k), T(k,k-1), etc.
% and O_k is the vector of Stevens operators O(k,q=k), O(k,k-1), ... O(k,0),...
% O(k,-k). Stevens operators with non-negative q are cosine tesseral operators,
% and those with negative q are sine tesseral operators.
%
% C is the (2k+1)x(2k+1) transformation matrix. Its only non-negative elements
% are on the diagonal and the antidiagonal.
%
% For the inverse transform, use inv(C): T_k = inv(C)*O_k

function C = isto2stev(k)

if ~isnumeric(k) || numel(k)~=1 || ~isreal(k) || k<0 || mod(k,1)
  error('Rank k must be a nonnegative integer (0, 1, 2, ...).');
end

% Generate ISTO vector and Stevens operator vector for rank k
% (ordering q = k, k-1, ..., -k)
J = k/2; % use smallest possible J (result C is independent of J)
for q = -k:k
  T{k+1-q} = isto(J,[k q]);
  O{k+1-q} = stev(J,[k,q]);
end

% Calculate elements of tranformation matrix via scalar products <O|T>
% (only matrix elements with qT=+-qO are non-zero)
C = zeros(2*k+1);
for t = 1:2*k+1
  traceTT = trace(T{t}'*T{t});
  for o = [t 2*k+2-t]
    C(o,t) = trace(O{o}.'*T{t})/traceTT;
  end
end
