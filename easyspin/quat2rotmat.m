%  quat2rotmat Convert a quaternion to a rotation matrix.
%
%  mat = quat2rotmat(q);
%
%  Input:
%     q          normalized quaternion
%
%  Output:
%     mat        rotation matrix

function mat = quat2rotmat(q)
qshape = size(q);
if numel(qshape)==2, qshape(3) = 1; end
    
if qshape(1) ~= 4
    error('Size of first dimension must equal 4.')
end
if any(1.0-sum(q.*q, 1) > 1e-5)
    error('Input is not normalized.')
end

q0 = q(1, :, :);
q1 = q(2, :, :);
q2 = q(3, :, :);
q3 = q(4, :, :);

mat = zeros(3, 3, qshape(2), qshape(3));

mat(1, 1, :, :) = 1-2*q2.*q2-2*q3.*q3;
mat(1, 2, :, :) = 2*q1.*q2-2*q3.*q0;
mat(1, 3, :, :) = 2*q1.*q3+2*q2.*q0;
mat(2, 1, :, :) = 2*q1.*q2+2*q3.*q0;
mat(2, 2, :, :) = 1-2*q1.*q1-2*q3.*q3;
mat(2, 3, :, :) = 2*q2.*q3-2*q1.*q0;
mat(3, 1, :, :) = 2*q1.*q3-2*q2.*q0;
mat(3, 2, :, :) = 2*q2.*q3+2*q1.*q0;
mat(3, 3, :, :) = 1-2*q1.*q1-2*q2.*q2;

% mat = 2*[0.5-q2.*q2-q3.*q3, q1.*q2-q3.*q0, q1.*q3+q2.*q0;
%          q1.*q2+q3.*q0, 0.5-q1.*q1-q3.*q3, q2.*q3-q1.*q0;
%          q1.*q3-q2.*q0, q2.*q3+q1.*q0, 0.5-q1.*q1-q2.*q2];

% mat = [q0^2+q1^2-q2^2-q3^2,     -2*(q0*q3-q1*q2),       2*(q0*q2+q1*q3);
%            2*(q0*q3+q1*q2),  q0^2-q1^2+q2^2-q3^2,      -2*(q0*q1-q2*q3);
%           -2*(q0*q2-q1*q3),      2*(q0*q1+q2*q3),   q0^2-q1^2-q2^2+q3^2];
end