% Calculates kron(eye(nI),B) for sparse B without multiplications
function K = speyekron(nI,B)

if ~issparse(B)
  error('B must be sparse.');
end

[mb,nb] = size(B);
ia = (1:nI).';
[ib,jb,sb] = find(B);
ib = ib(:); jb = jb(:); sb = sb(:);
ik = bsxfun(@plus, mb*(ia-1).', ib);
jk = bsxfun(@plus, nb*(ia-1).', jb);
K = sparse(ik,jk,repmat(sb,1,nI),nI*mb,nI*nb);

return
