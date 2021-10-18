%  quat2euler  Convert a unit quaternion into a sequence of Euler angles
%              using the z-y-z convention.
%
%   [alpha, beta, gamma] = quat2euler(q);
%
%   Input:
%     q              numeric, size = (4,...)
%                    normalized quaternion
%
%   Output:
%    alpha          double or numeric, size = (1,...)
%                   first Euler angle, in radians
%
%    beta           double or numeric, size = (1,...)
%                   second Euler angle, in radians
%
%    gamma          double or numeric, size = (1,...)
%                   third Euler angle, in radians

function varargout = quat2euler(varargin)

if (nargin==0), help(mfilename); return; end

switch nargin
  case 1
    % EasySpin uses passive transformations by default
%     q = varargin{1};
%     invert = 1;
    q = quatinv(varargin{1});
  case 2
    if ~ischar(varargin{2}) %|| ~isstring(varargin{2})
      error('If two arguments are given, the second must be a string.')
    end
    switch lower(varargin{2})
      case 'passive'
%         q = varargin{1};
%         invert = 1;
        q = quatinv(varargin{1});
      case 'active'
        q = varargin{1};
%         invert = 0;
      otherwise
        error('The second argument must specify ''active'' or ''passive''.')
    end
  otherwise
    error('Too many input arguments.')
end

qshape = size(q);
    
if qshape(1) ~= 4 || ~isnumeric(q)
  error('q must be an array of size (4,...)')
end

diff = abs(1.0-squeeze(sqrt(sum(q.*q, 1))));
if any(diff(:) > 1e-10)
  error('q is not normalized.')
end

Index = cell(1, ndims(q));
Index(:) = {':'};

idx = q(1,Index{2:end}) < 0;

if any(idx(:))
  q(:,idx) = -q(:,idx);
end
               
q0 = q(1,Index{2:end});
q1 = q(2,Index{2:end});
q2 = q(3,Index{2:end});
q3 = q(4,Index{2:end});

sy = 2*sqrt((q2.*q3+q1.*q0).^2 + (q1.*q3-q2.*q0).^2);

idx = sy < 1e-10;

alpha = atan2( 2.*(q2.*q3 ...
                 - q0.*q1), ...
               2.*(q1.*q3 ...
                 + q0.*q2) );
alpha(alpha<0) = alpha(alpha<0) + 2*pi;

% since beta can sometimes contain very small imaginary parts, take the
% real part
% beta = real( acos( q0.^2 - q1.^2 ...
%                  - q2.^2 + q3.^2) );
beta = real( atan2( sy, 1 - 2*q1.^2 - 2*q2.^2 ) );

gamma = atan2( 2.*(q2.*q3 ...
                 + q0.*q1), ...
              -2.*(q1.*q3 ...
                 - q0.*q2) );
gamma(gamma<0) = gamma(gamma<0) + 2*pi;
               
% if invert
%   temp = alpha;
%   alpha = -gamma;
%   beta = -beta;
%   gamma = -temp;
% end
               
% if any(idx)
%   alpha(idx) = zeros(size(alpha(idx)));
%   
%   beta(idx) = atan2( sy(idx), 1 - 2*q1(idx).^2 - 2*q2(idx).^2 );
%   
%   gamma(idx) = atan2( -2.*(-q1(idx).*q2(idx) ...
%                           + q3(idx).*q0(idx)), ...
%                    1 - 2.*( q1(idx).*q1(idx) ...
%                           - q3(idx).*q3(idx)) );
% end

switch (nargout)
  case 0
%     alpha, beta, gamma
    varargout = {cat(1, alpha, beta, gamma)};
  case 1
    varargout = {cat(1, alpha, beta, gamma)};
  case 3
    varargout = {alpha, beta, gamma};
  otherwise
    error('Incorrect number of output arguments.')
end

end