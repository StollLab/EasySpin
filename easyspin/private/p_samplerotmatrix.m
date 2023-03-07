% Parses Exp.SampleRotation, including letter codes and supplying a default if
% it's empty, and returns the (active) rotation matrix R.
%
%    Rrot = p_samplerotmatrix({rho,n})
%    vrot = Rrot*v
%
% vrot corresponds to the vector v (actively) rotated by angle rho around
% axis n in the counterclockwise direction.
%
% Simple example:
%   rotaxis = [0;0;1];
%   rho = 10*pi/180;
%   Rrot = rotaxi2mat(rotaxis,rho);
%   v = [1;0;0]
%   vrot = Rrot.'*v

function [R_samplerot,rotateSample] = p_samplerotmatrix(SampleRotation)

if isempty(SampleRotation)
  R_samplerot = eye(3);
  rotateSample = false;
  return
end

% Determine rotation axis and angle
if iscell(SampleRotation) && numel(SampleRotation)==2
  rho = SampleRotation{1};
  nRot = SampleRotation{2};
  if ischar(nRot)
    nRot = letter2vec(nRot);
  end
else
  error('Exp.SampleRotation must be of the form {rho,nL}, with the rotation angle rho and the lab-frame rotation axis nL.');
end
if numel(rho)~=1
  error('Exp.SampleRotation: the first element must a single rotation angle rho.');
end

% Calculate rotation matrix
R_samplerot = rotaxi2mat(nRot,rho).';

rotateSample = rho~=0;

end
