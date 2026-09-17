%  quatmult Perform quaternion multiplication.
%
%  t = quatmult(q, r);
%
%  Input:
%      q              numeric, size = (4,...)
%                     normalized quaternion
%
%      r              numeric, size = (4,...)
%                     normalized quaternion
%
%  Output:
%      t              numeric, size = (4,...)
%                     normalized quaternion

function t = quatmult(q, r)

qshape = size(q);
rshape = size(r);
    
if qshape(1) ~= 4 || rshape(1) ~= 4 || ~isnumeric(q) || ~isnumeric(r)
    error('q and r must be arrays of size (4,...)')
end

if ~isequal(size(q),size(r))
  error('Size input arrays must be equal.')
end

q_ = reshape(q,4,[]);
r_ = reshape(r,4,[]);

q0 = q_(1,:);
q1 = q_(2,:);
q2 = q_(3,:);
q3 = q_(4,:);

r0 = r_(1,:);
r1 = r_(2,:);
r2 = r_(3,:);
r3 = r_(4,:);

t_ = zeros(size(q_));

t_(1,:) = r0.*q0 - r1.*q1 - r2.*q2 - r3.*q3;
t_(2,:) = r0.*q1 + r1.*q0 - r2.*q3 + r3.*q2;
t_(3,:) = r0.*q2 + r1.*q3 + r2.*q0 - r3.*q1;
t_(4,:) = r0.*q3 - r1.*q2 + r2.*q1 + r3.*q0;

t = reshape(t_,qshape);

end
