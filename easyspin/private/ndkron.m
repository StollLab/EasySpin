%  ndkron  Kronecker product of two N-D arrays.
%
%          The size of dimensions 3,4,... need to be identical between A
%          and B.
%
%  C = ndkron(A,B);
%
%  Input:
%      A              numeric, size = (n,m,...)
%
%      B              numeric, size = (m,p,...)
%
%  Output:
%      C              numeric, size = (n,p,...)
%                     Kronecker product of A and B,
%                       $C = A \otimes B$

function C = ndkron(A,B)

if ~isnumeric(A) || ~isnumeric(B)
  error('Inputs must be of numeric type.')
end

if ndims(A)~=ndims(B)
  error('Inputs should have the same number of dimensions.')
end

if ismatrix(A) && ismatrix(B)
  % use built-in Matlab function if they are 2-D arrays
  C = kron(A,B);
  return
end

sizeA = size(A);
sizeB = size(B);

ma = sizeA(1);
na = sizeA(2);

mb = sizeB(1);
nb = sizeB(2);

if na~=mb
  error('size(A,2) and size(B,1) must be equal.')
elseif ~allclose(sizeA(3:end),sizeB(3:end))
  error('Input dimensions (3,4,...) must be equal.')
end

A = reshape(A,[1, ma, 1, na, sizeA(3:end)]);
B = reshape(B,[mb, 1, nb, 1, sizeB(3:end)]);
C = reshape(A.*B, [ma*mb, na*nb, sizeA(3:end)]);

end

