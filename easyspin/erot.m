% erot  Rotation/transformation matrix from Euler angles
%
%   Rp = erot(Angles)
%   Rp = erot(alpha,beta,gamma)
%   [xc,yc,zc] = erot(...,'cols')
%   [xr,yr,zr] = erot(...,'rows')
%
%   Computes a 3x3 rotation/transformation matrix Rp from
%   a vector of 3 Euler angles.
%
%   Input
%   - Angles: vector containing the three Euler angles (in radians) that define
%     the rotation. [alpha, beta, gamma] rotate the coordinate system around
%     [z,y',z''] counterclockwise in in that order.
%   - alpha, beta, gamma: the three Euler angles as defined above.
%   - 'cols', 'rows': tells erot to return either the three rows or the three
%     columns of the rotation matrix separately
%
%   Output
%   - Rp: matrix for the passive rotation/coordinate transformation
%       vec2 = Rp*vec1 or     % transform vector from frame 1 to frame 2
%       mat2 = Rp*mat1*Rp.'   % transform matrix from frame 1 to frame 2
%   - xc,yc,zc: three columns of the rotation matrix such that Rp = [xc yc zc]
%   - xr,yr,zr: three rows of the rotation matrix, as column vectors, such
%       that Rp = [xr.'; yr.'; zr.']

function varargout = erot(varargin)

if nargin==0, help(mfilename); return; end

switch nargin
  case {1,2}
    angles = varargin{1};
    if numel(angles)~=3
      error('  Three angles (either separately or in a 3-element array) expected.\n  You gave an %dx%d array with %d elements.',size(angles,1),size(angles,2),numel(angles));
    end
    gamma = angles(3);
    beta = angles(2);
    alpha = angles(1);
    if (nargin==2)
      option = varargin{2};
    else
      option = '';
    end
  case {3,4}
    gamma = varargin{3};
    beta = varargin{2};
    alpha = varargin{1};
    if (nargin==4)
      option = varargin{4};
    else
      option = '';
    end
  otherwise
    error('Wrong number of input arguments!');
end

if ~ischar(option)
    error('Last argument must be a string, either ''rows'' or ''cols''.')
end

switch option
  case ''
    returnRows = false;
  case 'rows'
    returnRows = true;
  case 'cols'
    returnRows = false;
  otherwise
    error('Last argument must be a string, either ''rows'' or ''cols''.')
end

if nargout==3 && isempty(option)
  error('Please specify whether erot() should return the 3 columns or the 3 rows of the matrix.');
end

if ~isempty(option) && nargout~=3
  error('3 outputs required if you specify ''rows'' or ''cols''.');
end

if ~any(nargout==[0 1 3])
  error('Wrong number of outputs!');
end

% Precalculate trigonometric functions of angles
sa = sin(alpha);
ca = cos(alpha);
sb = sin(beta);
cb = cos(beta);
sg = sin(gamma);
cg = cos(gamma);

% Compute passive rotation matrix
R = [ cg*cb*ca-sg*sa,   cg*cb*sa+sg*ca,  -cg*sb;...
     -sg*cb*ca-cg*sa,  -sg*cb*sa+cg*ca,   sg*sb;...
      sb*ca,            sb*sa,            cb];
%{
% explicitly:
Rg = [cg sg 0; -sg cg 0 ; 0 0 1];
Rb = [cb 0 -sb; 0 1 0 ; sb 0 cb];
Ra = [ca sa 0; -sa ca 0 ; 0 0 1];
R = Rg*Rb*Ra;
%}

if nargout==3
  if returnRows
    varargout = {R(1,:).',R(2,:).',R(3,:).'};
  else
    varargout = {R(:,1),R(:,2),R(:,3)};
  end
else
  varargout = {R};
end
