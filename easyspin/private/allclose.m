%  allclose  Check if ALL entries between two arrays agree within a 
%            specified tolerance.
%
%  val = allclose(A, B);
%  val = allclose(A, B, tol);
%
%   Input:
%     A              numeric
%
%     B              numeric
%
%
%   Output:
%     result         bool
%                    returns 1 if ALL entries in A-B are less than a 
%                    specified tolerance, 0 if not

% This function is inspired by "numpy.allclose" from the Numpy distribution.

function result = allclose(varargin)
% Check to see if two numeric arrays agree within a specified tolerance

% TODO add flexible functionality for NaNs

switch nargin
  case 2  % no tolerance specified, using built-in eps value
    A = varargin{1};
    B = varargin{2};
    tol = eps;
  case 3  % tolerance specified, use it
    A = varargin{1};
    B = varargin{2};
    tol = varargin{3};
    if ~isscalar(tol) || ~isa(tol,'double') || ~isreal(tol)
      error('Tol must be a real scalar of type double.')
    end
  otherwise
    error('Invalid number of inputs!')
end

if ~isnumeric(A) || ~isnumeric(B)
  error('Inputs A and B must be numeric arrays.')
end

diff = A - B;

NaNs = isnan(diff(:));
if any(NaNs), error('NaNs detected in difference!'); end

result = all(abs(diff(:)) < tol);

end

