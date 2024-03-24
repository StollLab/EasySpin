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
  R_L2S = cell(1,nSamples);
  for s = 1:nSamples
    R_L2S{s} = erot(Exp.SampleFrame(s,:));
  end
else
  R_L2S = {eye(3)};
  nSamples = 1;
end

% Check Exp.SampleRotation
activeSampleRotation = ~isempty(Exp.SampleRotation);
if activeSampleRotation
  if ~iscell(Exp.SampleRotation) || numel(Exp.SampleRotation)~=2
    error('Exp.SampleRotation must be of the form {nL,rho}, with the lab-frame rotation axis nL and the rotation angle rho.');
  end
end

% Apply (active) sample rotation if given
if activeSampleRotation
  % Get rotation axis and angle
  rotaxis = Exp.SampleRotation{1};
  rho = Exp.SampleRotation{2};
  if ischar(rotaxis)
    rotaxis = letter2vec(rotaxis);
  elseif numel(rotaxis)~=3
    error('Rotation axis (first element in Exp.SampleRotation) must be a 3-element vector or a letter code.');
  end
  nRotations = numel(rho);

  % Calculate rotation matrices
  % (vrot = Rrot*v corresponds to the vector v actively rotated by angle rho
  % around axis rotaxis in the counterclockwise direction.)
  for r = nRotations:-1:1
    Rrot{r} = rotaxi2mat(rotaxis,rho(r));
  end
  
  % Apply active sample rotations
  idx = 0;
  R_L2S0 = R_L2S;
  R_L2S = cell(1,nRotations*nSamples);
  for s = 1:nSamples
    for r = 1:nRotations
      idx = idx + 1;
      R_L2S{idx} = R_L2S0{s}*Rrot{r};
    end
  end

end 

% Convert to Euler angles
for idx = numel(R_L2S):-1:1
  angles_L2S(idx,:) = eulang(R_L2S{idx});
end

rotatedSample = any(angles_L2S(:));

end
