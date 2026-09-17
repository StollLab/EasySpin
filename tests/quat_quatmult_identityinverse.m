function ok = test()

% Check that quatmult() satisfies basic quaternion group properties:
%   q * qId  = q
%   qId * q  = q
%   q * qinv(q) = qId
% for a multi-dimensional array of quaternions.

N = 5;
M = 10;

q = rand(4,N,M);
q = q./sqrt(sum(q.*q,1));

qId = repmat([1;0;0;0],1,N,M);

diff1 = quatmult(q,qId) - q;
diff2 = quatmult(qId,q) - q;
diff3 = quatmult(q,quatinv(q)) - qId;

thr = 1e-10;
ok = all(abs(diff1(:))<thr) && all(abs(diff2(:))<thr) && all(abs(diff3(:))<thr);
