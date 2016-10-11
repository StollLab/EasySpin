%  quatmult Perform quaternion multiplication.
%
%  t = quatmult(q, r);
%
%  Input:
%     q          4x... array, normalized quaternion
%     r          4x... array, normalized quaternion
%
%  Output:
%     t          4x... array, normalized quaternion

function t = quatmult(q, r)
    
if or(size(q, 1) ~= 4, size(r, 1) ~= 4)
    error('Size of first dimension of inputs must equal 4.')
end

if ~isequal(size(q), size(r))
  error('Size input arrays must be equal.')
end

qlen = sum(q.*q, 1);
if any(1.0-qlen(:) > 1e-5)
    error('Input quaternion is not normalized.')
end

qIndex = cell(1, ndims(q));
qIndex(:) = {':'};

rIndex = cell(1, ndims(r));
rIndex(:) = {':'};

q0 = q(1,qIndex{2:end});
q1 = q(2,qIndex{2:end});
q2 = q(3,qIndex{2:end});
q3 = q(4,qIndex{2:end});

r0 = r(1,rIndex{2:end});
r1 = r(2,rIndex{2:end});
r2 = r(3,rIndex{2:end});
r3 = r(4,rIndex{2:end});


t0 = r0.*q0 - r1.*q1 - r2.*q2 - r3.*q3;
t1 = r0.*q1 + r1.*q0 - r2.*q3 + r3.*q2;
t2 = r0.*q2 + r1.*q3 + r2.*q0 - r3.*q1;
t3 = r0.*q3 - r1.*q2 + r2.*q1 + r3.*q0;

t = [t0; t1; t2; t3];

end