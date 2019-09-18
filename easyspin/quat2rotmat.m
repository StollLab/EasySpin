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
    
if qshape(1) ~= 4 || ~isnumeric(q)
  error('q must be an array of size (4,...)')
end

diff = 1-squeeze(sqrt(sum(q.*q, 1)));
if any(abs(diff(:)) > 1e-10)
%   plot(diff(1,:))
  error('Input is not normalized.')
end

Index = cell(1, ndims(q));
Index(:) = {':'};

q0 = q(1,Index{2:end});
q1 = q(2,Index{2:end});
q2 = q(3,Index{2:end});
q3 = q(4,Index{2:end});

% needed if 3rd dimension is singleton and 4th is not
mat = zeros([3,3,qshape(2:end)]); 

mat(1,1,Index{2:end}) = 1 - 2*q2.*q2 - 2*q3.*q3;
mat(1,2,Index{2:end}) = 2*q1.*q2 - 2*q0.*q3;
mat(1,3,Index{2:end}) = 2*q1.*q3 + 2*q0.*q2;
mat(2,1,Index{2:end}) = 2*q1.*q2 + 2*q0.*q3;
mat(2,2,Index{2:end}) = 1 - 2*q1.*q1 - 2*q3.*q3;
mat(2,3,Index{2:end}) = 2*q2.*q3 - 2*q0.*q1;
mat(3,1,Index{2:end}) = 2*q1.*q3 - 2*q0.*q2;
mat(3,2,Index{2:end}) = 2*q2.*q3 + 2*q0.*q1;
mat(3,3,Index{2:end}) = 1 - 2*q1.*q1 - 2*q2.*q2;

end