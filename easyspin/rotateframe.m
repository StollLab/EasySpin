% rotateframe   Rotate a frame around a given axis
%
%    ori_rot = rotateframe(Ori,nRot,rho)
%
%    Rotates a frame descibed by the three Euler angles in Ori around the axis
%    given in nRot by the rotation angles listed in rho.
%
%    Input:
%     Ori  Euler angles for the frame orientation
%     nRot   rotation axis
%               e.g., nRot = [1;0;0] is the x axis
%     rho     rotation angle, or list of rotation angles, for the rotation
%               around nRot
%
%    Output:
%     ori_rot  list of Euler angle sets for the rotated frames, one per row.
%
%    Example:
%       Ori0 = [0 45 0]*pi/180;
%       nRot = [1;0;0];
%       rho = (0:30:180)*pi/180;
%       Ori_list = rotateframe(nRot,rho);

function out = rotateframe(Frame,nRot,rho)

if nargin==0
  help(mfilename);
end

if numel(Frame)~=3
  error('First input (CryOri) must contain three numbers, the three Euler angles.');
end

if ischar(nRot)
  nRot = letter2vec(nRot);
else
  if numel(nRot)~=3
    error('Second input (nRot) must be either a 3-element vector or a letter indicating the vector (''x'', ''y'', etc).');
  end
end

xyzC_L = erot(Frame);
% xyzC_L col 1,2,3: crystal axis xC,yC,zC in lab frame coordinates
% nRot: rotation axis

skipFitting = true;

for irho = 1:numel(rho)
  R = rotaxi2mat(nRot,rho(irho));
  xyzC_L_rotated = R.'*xyzC_L;
  angles_rotated(irho,:) = eulang(xyzC_L_rotated,skipFitting);
end

out = angles_rotated;
