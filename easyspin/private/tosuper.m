%  tosuper  Form a superoperator.
%
%           Note: because of Matlab's column-major ordering of array
%           elements, vectorized operators using "A(:)" will be ordered by 
%           "column-stacking", i.e. 
%              [1 2; 3 4] => [1; 3; 2; 4],
%           so the order of Kronecker products here will be reversed with 
%           respect to the textbook expressions. For example, the 
%           Liouvillian, with I being the identity operator, would be 
%           written as follows:
%              here:     L = kron(H,I) - kron(I,H)  (coulmn-major ordering)
%              textbook: L = kron(I,H) - kron(H,I)  (row-major ordering)
%
%  Asuper = tosuper(A,form);
%
%  Input:
%      A              numeric, size = (n,n,...)
%                     lefthand-acting matrix
%
%      form           char
%                     indicate what kind of superoperator to generate,
%                     either "C" for commutator or "A" for anti-commutator
%
%  Output:
%      Asuper        numeric, size = (n^2,n^2,...)
%                     superoperator

function Asuper = tosuper(A,form)

if ~ischar(form)
  error('Form must be character array.')
end

form = upper(form);

if ~isnumeric(A)
  error('Input must of numeric type.')
end

sizeA = size(A);

rA = sizeA(1);
cA = sizeA(2);

if rA~=cA
  error('Sizes of first 2 dimensions must be equal.')
end

ndeye = repmat(eye(rA),[1,1,sizeA(3:end)]);

switch form
  case 'C'
    Asuper = ndkron(ndeye, A) ...
             - ndkron(permute(A,[2,1,3:ndims(A)]), ndeye);
  case 'A'
    Asuper = ndkron(ndeye, A) ...
             + ndkron(permute(A,[2,1,3:ndims(A)]), ndeye);
  otherwise
    error('Expected "C" or "A". Please check documentation.')
end

end

