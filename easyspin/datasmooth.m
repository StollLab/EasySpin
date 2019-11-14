% datasmooth  Moving average smoothing and differentiation 
%
%   yy = datasmooth(y,m);
%   yy = datasmooth(y,m,'binom')
%   yy = datasmooth(y,m,'flat')
%   yy = datasmooth(y,m,'savgol')
%   yy = datasmooth(y,m,'savgol',p)
%   yy = datasmooth(y,m,'savgol',p,dif)
%
%   Computes the 2*m+1 - point moving average of y.
%   If y is a matrix, datasmooth operates along columns.
%
%   yy(i) is a weighted average of y(i-m:i+m). y is
%   enlarged at start and end by its start and end values,
%   e.g. y(0) = y(-1) = y(1), y(end+2) = y(end).
%
%   If 'flat' is given, the moving average is unweighted.
%
%   If 'binom' is specified, binomial coefficients are
%   used as weighting factors. This is the default method.
%
%   If 'savgol' is specified, a least-squares smoothing using 
%   the Savitzky-Golay polynomial filter of order p is computed.
%   It least-square fits p-order polynomials to 2*m+1 wide frames.
%   If dif>0, a derivative of y is computed at the same time. E.g.
%   if dif=3, y is denoised and its third derivative is returned.
%
%   If the smoothing filter is omitted, 'binom' is the default.
%   If p is omitted, 2 is the default.
%   If dif is omitted, 0 is the default.

function y_filtered = datasmooth(y,m,Method,PolyOrder,Derivative)

if nargin==0, help(mfilename); return; end
if nargin<2, m = 3; end
if nargin<3, Method = 'binom'; end
if nargin<4, PolyOrder = 2; end
if nargin<5, Derivative = 0; end

if ~isreal(m) || numel(m)~=1 || ~isfinite(m) || (m<0) || mod(m,1)
  error('m (second parameter) must be a positive integer!');
end
if m==0, y_filtered = y; return; end
if PolyOrder<Derivative
  error('Polynomial order must not be smaller than the derivative index!');
end
if PolyOrder<1 || mod(PolyOrder,1) || numel(PolyOrder)~=1
  error('Polynomial order p must be a positive integer!');
end

% Symmetric frame
mLeft = m;
mRight = m;
n = mLeft + mRight + 1;

switch Method
  case 'flat'
    if (nargin<3) || (nargin>3), error('Wrong number of input arguments!'); end
    Weights = ones(1,n)/n;
  case 'binom'
    if (nargin<2) || (nargin>3), error('Wrong number of input arguments!'); end
    Weights = diag(flipud(pascal(n))).';
    Weights = Weights/2^(n-1);
  case 'savgol'
    X = repmat((-mLeft:mRight).',1,PolyOrder+1).^repmat(0:PolyOrder,n,1);
    F = pinv(X);
    Weights = (-1)^Derivative*F(Derivative+1,:);
  otherwise
    error('Unknown value for third argument!');
end

RowVec = isrow(y);
if RowVec
  y = y.';
end

% Enlarge vector(s) at beginning and end.
yend = y(end,:);
y_expanded = [y(ones(1,mLeft),:); y; yend(ones(1,mRight+1),:)];

% Apply filter.
y_filtered = filter(Weights,1,y_expanded);

% Chop to right size.
y_filtered = y_filtered(n:end-1,:);

if RowVec
  y_filtered = y_filtered.';
end

return
