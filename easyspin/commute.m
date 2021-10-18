% commute  Commutator of two matrices 
%
%   C = commute(A,B)
%
%   Returns the commutator [A,B]=A*B-B*A of the
%   two matrices A and B. A and B must be
%   square and of the same order.

function C = commute(A,B)

if (nargin==0), help(mfilename); return; end

if ~isnumeric(A) | ~isnumeric(B)
  error('A and B must be matrices!');
end

if any(size(A)~=size(B))
  error('A and B must have the same size!');
end

if size(A,1)~=size(A,2)
  error('A and B must be square matrices!');
end

C = A*B - B*A;

return
