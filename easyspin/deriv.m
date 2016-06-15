% deriv   Derivative of a vector or matrix 
%
%   dydx = deriv(y);
%   dydx = deriv(x,y);
%
%   Computes the derivate of the input vector.
%   If x is not given, x is assumed to be 1:length(y).
%   If y is a matrix, differentiation occurs along
%   columns.

function dydx = deriv(varargin)

switch nargin
case 0
  help(mfilename);
  return;
case 1
  y = varargin{1};
  x = [];
case 2
  x = varargin{1};
  y = varargin{2};
otherwise
  error('Wrong number of input arguments!');
end

RowVector = (numel(y)==size(y,2));
if (RowVector), y = y(:); end

if isempty(x)
  x = 1:size(y,1);
end

dydx = diff(y)./repmat(diff(x(:)),1,size(y,2));
dydx = (dydx([1 1:end],:)+dydx([1:end end],:))/2;

if (RowVector), dydx = dydx.'; end

return;
