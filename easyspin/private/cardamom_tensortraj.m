% given one or more trajectories of rotation matrices, calculate the 
% corresponding time dependence of an interaction tensor
%
% Input:
%   T          tensor, either 3-vector, or 3x3 matrix
%   RTraj      trajectories of rotation matrices, size (3,3,nTraj,nSteps)
%   RTrajInv   trajectories of rotation matrix inverses
%
% Output:
%   TTraj      tensor trajectories, size (3,3,nTraj,nSteps)

function TTraj = cardamom_tensortraj(T,RTraj,RTrajInv)

if isvector(T) && numel(T)==3
  T = diag(T);
elseif ismatrix(T) && all(size(T)==[3 3])
  % do nothing
else
  error('Interaction tensor T needs to be given as a 3-vector or 3x3 matrix.')
end

if size(RTraj,1)~=3 || size(RTraj,2)~=3
  error('Rotation matrix trajectory needs to have size 3x3xnTrajxnSteps.')
end

nTraj = size(RTraj,3);
nSteps = size(RTraj,4);

TTraj = multimatmult(RTraj, ...
                     multimatmult(repmat(T,1,1,nTraj,nSteps), RTrajInv));

end
