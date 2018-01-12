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

if ~isequal(size(q), size(r))
  error('Size input arrays must be equal.')
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

t = zeros(size(q));

t(1,qIndex{2:end}) = r0.*q0 - r1.*q1 - r2.*q2 - r3.*q3;
t(2,qIndex{2:end}) = r0.*q1 + r1.*q0 - r2.*q3 + r3.*q2;
t(3,qIndex{2:end}) = r0.*q2 + r1.*q3 + r2.*q0 - r3.*q1;
t(4,qIndex{2:end}) = r0.*q3 - r1.*q2 + r2.*q1 + r3.*q0;

% t = [t0; t1; t2; t3];

end