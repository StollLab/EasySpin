% rotatecrystal   Rotate a crystal around a given axis
%
%    ori_rot = rotatecrystal(CryOri,nRotL,rho)
%
%    Rotates the crystal frame descibed by the three Euler angles
%    in CryOri around the axis given in nRotL by the rotation angles
%    listed in rho.
%
%    Input:
%     CryOri  Euler angles for the crystal->lab frame transformation
%               (same as Exp.CrystalOrientation)
%     nRotL   rotation axis, in lab coordinates
%               e.g., nRotL = [1;0;0] is the lab x axis, xL
%     rho     rotation angle, or list of rotation angles, for the rotation
%               around nRotL
%
%    Output:
%     ori_rot  list of Euler angle sets for the rotated frames, one per row.
%                   (can be directly used as Exp.CrystalOrientation in simulations).
%
%    Example:
%       COri0 = [0 45 0]*pi/180;
%       nRotL = [1;0;0];
%       rho = (0:30:180)*pi/180;
%       COri_list = rotatecrystal(nRotL,rho);
%       Exp.CrystalOrientation = COri_list;

function out = rotatecrystal(angCrystalOrientation,nRotL,rho)

if nargin==0
  help(mfilename);
end

if numel(angCrystalOrientation)~=3
  error('First input (CryOri) must contain three numbers, the three Euler angles.');
end

if ischar(nRotL)
  nRotL = letter2vec(nRotL);
else
  if numel(nRotL)~=3
    error('Second input (nRot) must be either a 3-element vector or a letter indicating the vector (''x'', ''y'', etc).');
  end
end

xyzC_L = erot(angCrystalOrientation);
% xyzC_L col 1,2,3: crystal axis xC,yC,zC in lab frame coordinates
% nRot: rotation axis, in lab frame coordinates

skipFitting = true;

for irho = 1:numel(rho)
  Rot_L = rotaxi2mat(nRotL,rho(irho));
  xyzC_L_rotated = Rot_L.'*xyzC_L;
  angles_rotated(irho,:) = eulang(xyzC_L_rotated,skipFitting);
end

out = angles_rotated;
