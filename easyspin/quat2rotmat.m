%  quat2rotmat Convert a quaternion to a rotation matrix.
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

diff = abs(1.0-squeeze(sum(q.*q, 1)));
if any(diff(:) > 1e-13)
%   plot(diff(1,:))
  error('Input is not normalized.')
end

Index = cell(1, ndims(q));
Index(:) = {':'};

q0 = q(1,Index{2:end});
q1 = q(2,Index{2:end});
q2 = q(3,Index{2:end});
q3 = q(4,Index{2:end});

sq1 = q1.*q1;
sq2 = q2.*q2;
sq3 = q3.*q3;

q0q1 = q0.*q1;
q0q2 = q0.*q2;
q0q3 = q0.*q3;
q1q2 = q1.*q2;
q1q3 = q1.*q3;
q2q3 = q2.*q3;

% mat = zeros(3,3,qshape(2:end));

mat(1,1,Index{2:end}) = 1 - 2*sq2 - 2*sq3;
mat(1,2,Index{2:end}) = 2*q1q2 - 2*q0q3;
mat(1,3,Index{2:end}) = 2*q1q3 + 2*q0q2;
mat(2,1,Index{2:end}) = 2*q1q2 + 2*q0q3;
mat(2,2,Index{2:end}) = 1 - 2*sq1 - 2*sq3;
mat(2,3,Index{2:end}) = 2*q2q3 - 2*q0q1;
mat(3,1,Index{2:end}) = 2*q1q3 - 2*q0q2;
mat(3,2,Index{2:end}) = 2*q2q3 + 2*q0q1;
mat(3,3,Index{2:end}) = 1 - 2*sq1 - 2*sq2;

end