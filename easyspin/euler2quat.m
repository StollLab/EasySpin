%  euler2quat  Convert Euler angles to unit quaternions.
%
%  q = euler2quat(Angles);
%  q = euler2quat(alpha, beta, gamma);
%
%  Input:
%     Angles         numeric array, size = (3,...)
%                    the rows correspond to first, second and third Euler
%                    angle(s) (alpha, beta, and gamma, respectively), in radians
%     alpha          double or numeric array, size = (1,...)
%                    first Euler angle(s), in radians
%     beta           double or numeric array, size = (1,...)
%                    second Euler angle(s), in radians
%     gamma          double or numeric array, size = (1,...)
%                    third Euler angle(s), in radians
%
%  Output:
%     q              numeric, size = (4,...)
%                    normalized quaternion(s)

function q = euler2quat(varargin)

if nargin==0, help(mfilename); return; end

switch nargin
  case {1,2}  % 3-vector or 3x... shaped array
    angles = varargin{1};
    if isvector(angles)
      % 3-vector input
      if numel(angles)~=3
        error('Input must be 3 separate angles, a 3-vector, or a (3,...) shaped array.')
      end
      alpha = angles(1);
      beta = angles(2);
      gamma = angles(3);
    elseif size(angles,1)==3
      % 3x... shaped array input
      idx = cell(1, ndims(angles));
      idx(:) = {':'};
      alpha = angles(1,idx{2:end});
      beta = angles(2,idx{2:end});
      gamma = angles(3,idx{2:end});
    else
      error('Input must be 3 separate angles, a 3-vector, or a (3,...) shaped array.')
    end
    if (nargin==2)
      option = varargin{2};
    else
      option = '';
    end
  case {3,4} % 3 separate angles or 3 1-D arrays
    if ~isvector(varargin{1}) || ~isvector(varargin{2}) || ~isvector(varargin{1})
      error('If three arguments are given, they must be scalars or vectors.')
    end
    alpha = varargin{1};
    beta = varargin{2};
    gamma = varargin{3};
    if iscolumn(alpha), alpha = alpha.'; end
    if iscolumn(beta), beta = beta.'; end
    if iscolumn(gamma), gamma = gamma.'; end
    if (nargin==4)
      option = varargin{4};
    else
      option = '';
    end

end

if numel(alpha) ~= numel(beta) || numel(beta) ~= numel(gamma)
  error('Inputs must have the same number of elements.')
end

switch lower(option)
  case ''
    % EasySpin uses passive transformations by default
    temp = alpha;
    alpha = -gamma;
    beta = -beta;
    gamma = -temp;
  case 'passive'
    temp = alpha;
    alpha = -gamma;
    beta = -beta;
    gamma = -temp;
  case 'active'
    % do nothing
  otherwise
    error('Last argument must be a string, either ''passive'' or ''active''.')
end

q0 = cos(beta/2).*cos((gamma+alpha)/2);
q1 = sin(beta/2).*sin((gamma-alpha)/2);
q2 = sin(beta/2).*cos((gamma-alpha)/2);
q3 = cos(beta/2).*sin((gamma+alpha)/2);

q = [q0; q1; q2; q3];

% since q and -q give the same transformation, arbitrarily define q0 as 
% always positive

Index = cell(1, ndims(q));
Index(:) = {':'};

idx = q(1,Index{2:end}) < 0;

if any(idx(:))
  q(:,idx) = -q(:,idx);
end

end