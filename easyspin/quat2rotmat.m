%  quat2rotmat Convert unit quaternions to rotation matrices.
%
%  mat = quat2rotmat(q);
%
%  Input:
%      q              numeric, size = (4,...)
%                     normalized quaternion
%
%  Output:
%      mat            numeric, size = (3,3,...)
%                     rotation matrix

function mat = quat2rotmat(q)

qshape = size(q);
    
if ~isnumeric(q) || qshape(1)~=4
  error('q must be an array of size (4,...)')
end

diff = 1 - squeeze(sqrt(sum(q.*q, 1)));
if any(abs(diff(:)) > 1e-10)
  error('Input is not normalized.')
end

idx = repmat({':'},1,ndims(q)-1);

q0 = q(1,idx{:});
q1 = q(2,idx{:});
q2 = q(3,idx{:});
q3 = q(4,idx{:});

% needed if 3rd dimension is singleton and 4th is not
mat = zeros([3,3,qshape(2:end)]); 

mat(1,1,idx{:}) = 1 - 2*q2.*q2 - 2*q3.*q3;
mat(1,2,idx{:}) = 2*q1.*q2 - 2*q0.*q3;
mat(1,3,idx{:}) = 2*q1.*q3 + 2*q0.*q2;
mat(2,1,idx{:}) = 2*q1.*q2 + 2*q0.*q3;
mat(2,2,idx{:}) = 1 - 2*q1.*q1 - 2*q3.*q3;
mat(2,3,idx{:}) = 2*q2.*q3 - 2*q0.*q1;
mat(3,1,idx{:}) = 2*q1.*q3 - 2*q0.*q2;
mat(3,2,idx{:}) = 2*q2.*q3 + 2*q0.*q1;
mat(3,3,idx{:}) = 1 - 2*q1.*q1 - 2*q2.*q2;

end
