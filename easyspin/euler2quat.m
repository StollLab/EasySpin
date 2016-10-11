%  euler2quat  Convert set of Euler angles to a quaternion.
%
%  q = euler2quat(alpha, beta, gamma);
%
%  Input:
%     alpha          double or 1xN array, first Euler angle
%     beta           double or 1xN array, second Euler angle
%     gamma          double or 1xN array, third Euler angle
%
%  Output:
%     q              4xN array, normalized quaternion

function q = euler2quat(varargin)

if (nargin==0), help(mfilename); return; end

switch (nargin)
  case {1,2}  % 3-vector or 3x... shaped array
    angles = varargin{1};
    if isvector(angles)
      % 3-vector input
      if numel(angles)~=3
        error('Input must be 3 separate angles, a 3-vector, or a 3x... shaped array.')
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
      error('Input must be 3 separate angles, a 3-vector, or a 3x... shaped array.')
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
    if size(alpha,1)>1, alpha = alpha.'; end
    if size(beta,1)>1, beta = beta.'; end
    if size(gamma,1)>1, gamma = gamma.'; end
    if (nargin==4)
      option = varargin{4};
    else
      option = '';
    end

if numel(alpha) ~= numel(beta) || numel(beta) ~= numel(gamma)
  error('Inputs must have the same number of elements.')
end

switch option
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
    % Do nothing
  otherwise
    error('Last argument must be a string, either ''passive'' or ''active''.')
end

q0 = cos(beta/2).*cos((gamma+alpha)/2);
q1 = sin(beta/2).*sin((gamma-alpha)/2);
q2 = sin(beta/2).*cos((gamma-alpha)/2);
q3 = cos(beta/2).*sin((gamma+alpha)/2);

q = [q0; q1; q2; q3];

end