% given one or more trajectories of rotation matrices, calculate the 
% corresponding time dependence of an interaction tensor

function TTraj = cardamom_tensortraj(T,RTraj,RTrajInv)

if numel(T)==3
  T = diag(T);
elseif size(T,1)==3&&size(T,2)==3
  % do nothing
else
  error('Interaction tensor T needs to be given as a 3-vector or 3x3 matrix.')
end

nTraj = size(RTraj,3);
nSteps = size(RTraj,4);

TTraj = multimatmult(RTraj, ...
                     multimatmult(repmat(T,1,1,nTraj,nSteps), RTrajInv));

end