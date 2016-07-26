%  euler2rotmat  Convert set of Euler angles to a rotation matrix.
%
%  mat = euler2rotmat(alpha, beta, gamma);
%
%  Input:
%     alpha          first Euler angle
%     beta           second Euler angle
%     gamma          third Euler angle
%
%  Output:
%     mat            rotation matrix

function mat = euler2rotmat(alpha, beta, gamma)

q = euler2quat(alpha, beta, gamma);

mat = quat2rotmat(q);

end