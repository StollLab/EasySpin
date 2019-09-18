%  rotmat2quat Convert rotation matrices to unit quaternion.
%
%  q = rotmat2quat(mat);
%
%  Input:
%      mat            numeric, size = (3,3,...)
%                     rotation matrix
%
%  Output:
%      q              numeric, size = (4,...)
%                     normalized quaternion

% Implementation based on:
%   http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm

function q = rotmat2quat(mat)
    
if size(mat,1)~=3||size(mat,2)~=3||~isnumeric(mat)
  error('mat must be an array of size (3,3,...)')
end

if ndims(mat)>4
  error('Arrays with greater than 4 dimensions are not supported.')
end

if ismatrix(mat)
  diff = abs(eye(3)-mat*mat');
else
  matInv = permute(mat,[2,1,3,4]);
  diff = abs(repmat(eye(3),1,1,size(mat,3),size(mat,4))...
             - multimatmult(mat,matInv));
end

if any(diff(:) > 1e-10)
  error('The rotation matrices are not orthogonal.')
end

% Index = cell(1, ndims(mat));
% Index(:) = {':'};
q = zeros(4,size(mat,3),size(mat,4));

tr = mat(1,1,:,:) + mat(2,2,:,:) + mat(3,3,:,:);

for k=1:size(mat,4)
  for n=1:size(mat,3)
    imat = mat(:,:,n,k);
    itr = tr(:,:,n,k);
    m11 = imat(1,1);
    m12 = imat(1,2);
    m13 = imat(1,3);
    m21 = imat(2,1);
    m22 = imat(2,2);
    m23 = imat(2,3);
    m31 = imat(3,1);
    m32 = imat(3,2);
    m33 = imat(3,3);
    if itr>0
      S = 2*sqrt(itr+1);  % S = 4*qw
      qw = S/4;
      qx = (m32 - m23)/S;
      qy = (m13 - m31)/S;
      qz = (m21 - m12)/S;
    elseif (m11>m22)&&(m11>m33)
      S = 2*sqrt(1 + m11 - m22 - m33);  % S = 4*qx
      qw = (m32 - m23)/S;
      qx = S/4;
      qy = (m12 + m21)/S;
      qz = (m13 + m31)/S;
    elseif m22>m33
      S = 2*sqrt(1 + m22 - m11 - m33);  % S = 4*qy
      qw = (m13 - m31)/S;
      qx = (m12 + m21)/S;
      qy = S/4;
      qz = (m23 + m32)/S;
    else
      S = 2*sqrt(1 + m33 - m11 - m22);  % S = 4*qz
      qw = (m21 - m12)/S;
      qx = (m13 + m31)/S;
      qy = (m23 + m32)/S;
      qz = S/4;
    end
    q(:,n,k) = [qw; qx; qy; qz];
  end
end

diff = abs(1.0-squeeze(sum(q.*q, 1)));
if any(diff(:) > 1e-10)
%   plot(diff(1,:))
  error('Output q is not normalized.')
end

end