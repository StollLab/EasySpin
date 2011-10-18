% erot  Rotation matrix from Euler angles
%
%   Rp = erot(Angles)
%   Rp = erot(alpha,beta,gamma)
%
%   Computes a 3x3 rotation matrix Rp from
%   a vector of 3 Euler angles.
%
%   Input
%   - Angles: vector containing the three
%     Euler angles (in radians) that define
%     the rotation. [alpha, beta, gamma] rotate
%     the coordinate system around [z,y',z'']
%     counterclockwise in in that order.
%   - alpha, beta, gamma: the three Euler angles
%     as defined above.
%
%   Output
%   - Rp: matrix for the passive rotation
%       vec1 = Rp*vec or
%       mat1 = Rp*mat*Rp.'

% Syntax with 3 *output* arguments is undocumented.
% Don't remove, is used by resfields and endorfrq.

function varargout = erot(varargin)

if (nargin==0), help(mfilename); return; end

switch (nargin)
  case 1
    angles = varargin{1};
    if numel(angles)~=3
      error('Three angles expected.');
    end
    gamma = angles(3);
    beta = angles(2);
    alpha = angles(1);
  case 3
    gamma = varargin{3};
    beta = varargin{2};
    alpha = varargin{1};
  otherwise
    error('Wrong number of input arguments!');
end

switch (nargout)
  case 3
  case 1
  case 0
  otherwise
    error('Wrong number of outputs!');
end

% Precalculate trigonometric functions of angles
sa = sin(alpha); ca = cos(alpha);
sb = sin(beta);  cb = cos(beta);
sg = sin(gamma); cg = cos(gamma);

% Compute passive rotation matrix
R = [ cg*cb*ca-sg*sa,   cg*cb*sa+sg*ca,  -cg*sb;...
  -sg*cb*ca-cg*sa,  -sg*cb*sa+cg*ca,   sg*sb;...
  sb*ca,            sb*sa,            cb];

if (nargout==3)
  varargout = {R(1,:),R(2,:),R(3,:)};
else
  varargout = {R};
end
