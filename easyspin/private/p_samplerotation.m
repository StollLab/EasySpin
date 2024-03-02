% Calcualtes transformation matrix from lab frame to sample/crystal frame
% from Exp.SampleFrame and Exp.SampleRotation.

function [R_L2S,rotatedSample] = p_samplerotation(Exp)

% Determine lab-to-sample/crystal frame transformation(s), R_L2S
% - R_L2S col 1,2,3: lab axis 1,2,3 represented in sample/crystal frame
% - R_L2S row 1,2,3: crystal axis 1,2,3 represented in lab frame
if ~isempty(Exp.SampleFrame)
  nSamples = size(Exp.SampleFrame,1);
  logmsg(1,'  %d sample orientation(s) given in Exp.SampleFrame',nSamples);

  % Construct transformation matrices (lab frame to sample/crystal frame)
  for s = nSamples:-1:1
    R_L2S{s} = erot(Exp.SampleFrame(s,:));
  end
else
  R_L2S = {eye(3)};
  nSamples = 1;
end

% Apply (active) sample rotation if given
activeSampleRotation = isfield(Exp,'SampleRotation') && ~isempty(Exp.SampleRotation);
if activeSampleRotation
  Rrot = samplerotmatrix(Exp.SampleRotation);
  for s = nSamples:-1:1
    R_L2S{s} = R_L2S{s}*Rrot;
  end
end

% Convert to Euler angles
for s = nSamples:-1:1
  angles_L2S(s,:) = eulang(R_L2S{s});
end

rotatedSample = any(angles_L2S(:));

end


% Parses Exp.SampleRotation, including letter codes and supplying a default if
% it's empty, and returns the (active) rotation matrix R.
%
%    Rrot = samplerotmatrix({rho,n})
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


function R_samplerot = samplerotmatrix(SampleRotation)

if isempty(SampleRotation)
  R_samplerot = eye(3);
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
R_samplerot = rotaxi2mat(nRot,rho);

end
