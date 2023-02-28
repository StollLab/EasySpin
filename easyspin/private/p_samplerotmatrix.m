% Parses Exp.SampleRotation, including letter codes and supplying a default if
% it's empty, and returns the (active) rotation matrix.

function [R_sample,rotateSample] = p_samplerotmatrix(SampleRotation)

if isempty(SampleRotation)
  R_sample = eye(3);
  rotateSample = false;
  return
end

% Determine rotation axis and angle
if isnumeric(SampleRotation)
  rho = SampleRotation;
  nRot = [1;0;0]; % lab x axis (xL_L) as default
elseif iscell(SampleRotation) && numel(SampleRotation)==2
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
R_sample = rotaxi2mat(nRot,rho);

rotateSample = rho~=0;

end
