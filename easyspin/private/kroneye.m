% kroneye Kronecker product with identity matrix
%
%   K = kroneye(A)
%   K = kroneye(A,nI)
%
%   Calculates kron(A,eye(nI)).
%
%   If nI is not given, it is set equal to the dimension of A, which must be
%   square in that case.
%
%   The function works for both full and sparse matrices.

function K = kroneye(varargin)

manualMethod = false;

switch nargin
  case 1
    A  = varargin{1};
    nI = size(A,1);
    if size(A,2)~=nI
      error('Input matrix must be square.');
    end
  case 2
    A = varargin{1};
    nI = varargin{2};
  case 3
    A = varargin{1};
    nI = varargin{2};
    manualMethod = varargin{3};
end

if ~isscalar(nI)
  error('Second input (nI) must be square');
end

if manualMethod  
  K = manualkroneye(A,nI);  
else
  if issparse(A), I = speye(nI); else, I = eye(nI); end
  K = kron(A,I);
end

function K = manualkroneye(A,nI)

if issparse(A)
  
  [ma,na] = size(A);
  [ia,ja,sa] = find(A);
  ib = (1:nI).';
  ia = ia(:); ja = ja(:); sa = sa(:);
  ik = bsxfun(@plus, nI*(ia-1).', ib);
  jk = bsxfun(@plus, nI*(ja-1).', ib);
  K = sparse(ik,jk,repmat(sa.',nI,1),ma*nI,na*nI);
  
else
  
  [r,c] = size(A);
  K = zeros(r*nI,c*nI);
  Rows = 0:nI:(r-1)*nI;
  Cols = 0:nI:(c-1)*nI;
  for k = 1:nI
    Rows = Rows + 1;
    Cols = Cols + 1;
    K(Rows,Cols) = A;
  end
  
end
