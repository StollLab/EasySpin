% Calculates kron(A,eye(nI)) for sparse A without multiplications
function K = spkroneye(A,nI)

if ~issparse(A)
  error('A must be sparse.');
end

[ma,na] = size(A);
[ia,ja,sa] = find(A);
ib = (1:nI).';
ia = ia(:); ja = ja(:); sa = sa(:);
ik = bsxfun(@plus, nI*(ia-1).', ib);
jk = bsxfun(@plus, nI*(ja-1).', ib);
K = sparse(ik,jk,repmat(sa.',nI,1),ma*nI,na*nI);

return
