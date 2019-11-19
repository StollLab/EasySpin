% eyekron Kronecker product with identity matrix
%
%   K = eyekron(A)
%   K = eyekron(nI,A);
%
%   Calculates kron(eye(nI),A) without multiplications. This is signficantly
%   faster for large matrices.
%
%   If nI is not given, it is set equal to the dimension of A, which must be
%   square in that case.
%
%   The function works for both full and sparse matrices.

function K = eyekron(varargin)

switch nargin
  case 1
    A  = varargin{1};
    nI = size(A,1);
    if size(A,2)~=nI
      error('Input matrix must be square.');
    end
  case 2
    nI = varargin{1};
    A = varargin{2};
end

if ~issparse(A)
  
  [mb,nb] = size(A);
  ia = (1:nI).';
  [ib,jb,sb] = find(A);
  ib = ib(:); jb = jb(:); sb = sb(:);
  ik = bsxfun(@plus, mb*(ia-1).', ib);
  jk = bsxfun(@plus, nb*(ia-1).', jb);
  K = sparse(ik,jk,repmat(sb,1,nI),nI*mb,nI*nb);
  
else
  
  [r,c] = size(A);
  K = zeros(r*nI,c*nI);
  Rows = 1:r;
  Cols = 1:c;
  for k = 1:nI
    K(Rows,Cols) = A;
    Rows = Rows + r;
    Cols = Cols + c;
  end
  
end


return
