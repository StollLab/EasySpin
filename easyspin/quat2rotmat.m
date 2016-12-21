%  quat2rotmat Convert a quaternion to a rotation matrix.
%
%  mat = quat2rotmat(q);
%
%  Input:
%     q          4x... array, normalized quaternion
%
%  Output:
%     mat        3x3x... array, rotation matrix

function mat = quat2rotmat(q)
qshape = size(q);
    
if qshape(1) ~= 4 || ~isnumeric(q)
  error('q must be an array of size 4x...')
end

diff = abs(1.0-squeeze(sum(q.*q, 1)));
if any(diff(:) > 1e-13)
  plot(diff(1,:))
  error('Input is not normalized.')
end

Index = cell(1, ndims(q));
Index(:) = {':'};

q0 = q(1,Index{2:end});
q1 = q(2,Index{2:end});
q2 = q(3,Index{2:end});
q3 = q(4,Index{2:end});

% mat = zeros(3,3,Index{2:end});

mat(1,1,Index{2:end}) = 1 - 2*q2.*q2 - 2*q3.*q3;
mat(1,2,Index{2:end}) = 2*q1.*q2 - 2*q3.*q0;
mat(1,3,Index{2:end}) = 2*q1.*q3 + 2*q2.*q0;
mat(2,1,Index{2:end}) = 2*q1.*q2 + 2*q3.*q0;
mat(2,2,Index{2:end}) = 1 - 2*q1.*q1 - 2*q3.*q3;
mat(2,3,Index{2:end}) = 2*q2.*q3 - 2*q1.*q0;
mat(3,1,Index{2:end}) = 2*q1.*q3 - 2*q2.*q0;
mat(3,2,Index{2:end}) = 2*q2.*q3 + 2*q1.*q0;
mat(3,3,Index{2:end}) = 1 - 2*q1.*q1 - 2*q2.*q2;

end