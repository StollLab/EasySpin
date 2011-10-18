function [p,S,mu] = polyfitv(x,y,n)
% polyfitv    Fit polynomial to data, vectorized
%
%   p = polyfit(x,y,n);
%
%   polyfit finds the coefficients of a polynomial p(x) of
%   degree n that fits the data y best in a least-squares sense. p is a
%   row vector of length n+1 containing the polynomial coefficients in
%   descending powers, p(1)*x^n + p(2)*x^(n-1) +...+ p(n)*x + p(n+1).
%
%   If y is a matrix, polyfitv() fits a polynomial to each column. p
%   is an array with one polynomial in each row.

% The regression problem is formulated in matrix format as:
%
%    y = V*p    or
%
%          3  2
%    y = [x  x  x  1] [p3
%                      p2
%                      p1
%                      p0]
%
% where the vector p contains the coefficients to be found.  For a
% 7th order polynomial, matrix V would be:
%
% V = [x.^7 x.^6 x.^5 x.^4 x.^3 x.^2 x ones(size(x))];

x = x(:);
if numel(y)==length(y), y = y(:); end

if (nargout>2)
  mu = [mean(x); std(x)];
  x = (x - mu(1))/mu(2);
end

if size(y,1)~=numel(x)
  error('y has wrong size!');
end

% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1);
for j = n:-1:1
  V(:,j) = x.*V(:,j+1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);
ws = warning('off'); 
p = R\(Q'*y);    % Same as p = V\y;
warning(ws);
if size(R,2) > size(R,1)
   warning('MATLAB:polyfit:PolyNotUnique', ...
       'Polynomial is not unique; degree >= number of data points.')
elseif condest(R) > 1.0e10
    if nargout > 2
        warning('MATLAB:polyfit:RepeatedPoints', ...
            'Polynomial is badly conditioned. Remove repeated data points.')
    else
        warning('MATLAB:polyfit:RepeatedPointsOrRescale', ...
            ['Polynomial is badly conditioned. Remove repeated data points\n' ...
            '         or try centering and scaling as described in HELP POLYFIT.'])
    end
end
r = y - V*p;
p = p.';          % Polynomial coefficients are row vectors by convention.

% S is a structure containing three elements: the triangular factor from a
% QR decomposition of the Vandermonde matrix, the degrees of freedom and
% the norm of the residuals.
S.R = R;
S.df = length(x) - (n+1);
S.normr = norm(r);
