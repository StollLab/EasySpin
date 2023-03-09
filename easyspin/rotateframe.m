% rotateframe   Rotate a frame around a given axis
%
%    ori_rot = rotateframe(ori,nRot,rho)
%
%    Rotates a frame described by the three Euler angles in ori around the axis
%    given in nRot by the rotation angles listed in rho.
%
%    Input:
%     ori    Euler angles describing the frame orientation
%     nRot   rotation axis
%               e.g., nRot = [1;0;0] is the x axis
%               does not need to be normalized 
%     rho    rotation angle, or list of rotation angles, for the rotation
%               around nRot
%
%    Output:
%     ori_rot  list of Euler angle sets for the rotated frames, one per row.
%
%    Example:
%       ori0 = [0 45 0]*pi/180;
%       nRot = [1;0;0];
%       rho = (0:30:180)*pi/180;
%       ori_rot = rotateframe(ori0,nRot,rho);

function ori_rot = rotateframe(ori,nRot,rho)

if nargin==0
  help(mfilename);
end

if numel(ori)~=3
  error('First input (initial frame) must contain three numbers, the three Euler angles.');
end

if ischar(nRot)
  nRot = letter2vec(nRot);
else
  if numel(nRot)~=3
    error('Second input (nRot) must be either a 3-element vector or a letter indicating the vector (''x'', ''y'', etc).');
  end
end

% Get frame unit vectors
xyz = erot(ori).';  % transpose!
% columns of xyz are the x, y and z vector in the reference frame

skipFitting = true;

nrho = numel(rho);
ori_rot = zeros(nrho,3);
for irho = 1:numel(rho)
  R = rotaxi2mat(nRot,rho(irho));
  xyz_rotated = R.'*xyz;  % active rotation!
  ori_rot(irho,:) = eulang(xyz_rotated.',skipFitting);
end

end
