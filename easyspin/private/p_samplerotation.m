% Calcualtes transformation matrix from lab frame to sample/crystal frame
% from Exp.SampleFrame and Exp.SampleRotation.
%
%  Exp.SampleFrame    array of three Euler angles defining lab->sample frame
%                     transformation
%  Exp.SampleRotation {nRot,rho} active rotation of sample frame around
%                     axis nRot by angle(s) rho, in radians
%
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
  idx = 0;
  for s = nSamples:-1:1
    for r = numel(Rrot):-1:1
      idx = idx + 1;
      R_L2S{idx} = R_L2S{s}*Rrot{r};
    end
  end
end

% Convert to Euler angles
for s = numel(R_L2S):-1:1
  angles_L2S(s,:) = eulang(R_L2S{s});
end

rotatedSample = any(angles_L2S(:));

end


% Parses Exp.SampleRotation, including letter codes and supplying a default if
% it's empty, and returns the (active) rotation matrix R.
%
%    Rrot = samplerotmatrix({nRot,rho})
%    vrot = Rrot*v
%
% vrot corresponds to the vector v (actively) rotated by angle rho around
% axis nRot in the counterclockwise direction.
%
% Simple example:
%   nRot = [0;0;1];
%   rho = 10*pi/180;
%   Rrot = rotaxi2mat(nRot,rho);
%   v = [1;0;0]
%   vrot = Rrot.'*v


function R_samplerot = samplerotmatrix(SampleRotation)

if isempty(SampleRotation)
  R_samplerot = eye(3);
  return
end

% Determine rotation axis and angle
if iscell(SampleRotation) && numel(SampleRotation)==2
  nRot = SampleRotation{1};
  rho = SampleRotation{2};
  if ischar(nRot)
    nRot = letter2vec(nRot);
  elseif numel(nRot)~=3
    error('Rotation axis (first element in Exp.SampleRotation) must be a 3-element vector or a letter code.');
  end
else
  error('Exp.SampleRotation must be of the form {nL,rho}, with the lab-frame rotation axis nL and the rotation angle rho.');
end

% Calculate rotation matrix
for i = numel(rho):-1:1
  R_samplerot{i} = rotaxi2mat(nRot,rho(i));
end

end
