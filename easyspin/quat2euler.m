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
%                   first Euler angle
%
%    beta           double or numeric, size = (1,...)
%                   second Euler angle
%
%    gamma          double or numeric, size = (1,...)
%                   third Euler angle

%   References
%   ----------
%   [1] Mathworks Aerospace Toolbox

function [alpha, beta, gamma] = quat2euler(varargin)

if (nargin==0), help(mfilename); return; end

switch nargin
  case 1
    % EasySpin uses passive transformations by default
    q = quatinv(varargin{1});
  case 2
    if ~ischar(varargin{2}) %|| ~isstring(varargin{2})
      error('If two arguments are given, the second must be a string.')
    end
    switch lower(varargin{2})
      case 'passive'
        q = quatinv(varargin{1});
      case 'active'
        q = varargin{1};
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

diff = abs(1.0-squeeze(sum(q.*q, 1)));
if any(diff(:) > 1e-13)
  error('q is not normalized.')
end

Index = cell(1, ndims(q));
Index(:) = {':'};

alpha = squeeze( atan2( 2.*(q(3,Index{2:end}).*q(4,Index{2:end}) - q(1,Index{2:end}).*q(2,Index{2:end})), ...
                        2.*(q(2,Index{2:end}).*q(4,Index{2:end}) + q(1,Index{2:end}).*q(3,Index{2:end})))...
               );

beta = squeeze( acos( q(1,Index{2:end}).^2 - q(2,Index{2:end}).^2 - q(3,Index{2:end}).^2 + q(4,Index{2:end}).^2)...
              );

gamma = squeeze( atan2( 2.*(q(3,Index{2:end}).*q(4,Index{2:end}) + q(1,Index{2:end}).*q(2,Index{2:end})), ...
                       -2.*(q(2,Index{2:end}).*q(4,Index{2:end}) - q(1,Index{2:end}).*q(3,Index{2:end})))...
               );

end