% rotateframe   Rotate a frame around a given axis
%
%    ang = rotateframe(ang0,nRot,rho)
%    [ang,R] = rotateframe(ang0,nRot,rho)
%
%    Rotates a frame described by the three Euler angles in ang0 around the axis
%    given in nRot by the rotation angles listed in rho.
%
%    Input:
%     ang0   Euler angles describing the frame orientation, in radians
%     nRot   letter or vector specifying the rotation axis
%              nRot = [1;0;0]
%              nRot = 'x' is the x axis
%              vectors do not need to be normalized 
%     rho    rotation angle, or array of rotation angles, for the rotation
%               around nRot, in radians
%
%    Output:
%     ang    list of Euler angle sets for the rotated frames, one per row.
%     R      list of rotation matrices for the rotated frames, one per row.
%
%    Example:
%       ang0 = [0 45 0]*pi/180;
%       nRot = [1;0;0];
%       rho = (0:30:180)*pi/180;
%       ang = rotateframe(ang0,nRot,rho);

function [ang,R] = rotateframe(ang0,nRot,rho)

if nargin==0
  help(mfilename);
end

if numel(ang0)~=3
  error('First input (initial frame) must contain three numbers, the three Euler angles.');
end

if ischar(nRot)
  nRot = letter2vec(nRot);
else
  if numel(nRot)~=3
    error('Second input (nRot) must be either a 3-element vector or a letter indicating the vector (''x'', ''y'', etc).');
  end
end

% Preallocate outputs
nrho = numel(rho);
ang = zeros(nrho,3);
R = cell(nrho,1);

R0 = erot(ang0);
for irho = 1:numel(rho)
  Rrot = rotaxi2mat(nRot,rho(irho));  % matrix for active rotation
  R{irho} = R0*Rrot;
  ang(irho,:) = eulang(R{irho},true);
end

if isscalar(R)
  R = R{1};
end

end
