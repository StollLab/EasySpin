function ok = test()

% Check basic properties of quatinv():
%   - quatinv(q) is the quaternion conjugate (negates the vector part)
%   - quatinv(quatinv(q)) == q
%   - quatmult(q,quatinv(q)) == identity quaternion, for normalized q
% for a multi-dimensional array of quaternions.

N = 5;
M = 10;

q = rand(4,N,M);
q = q./sqrt(sum(q.*q,1));

qinv = quatinv(q);

diffConj = qinv - [q(1,:,:); -q(2,:,:); -q(3,:,:); -q(4,:,:)];

diffDouble = quatinv(qinv) - q;

qId = repmat([1;0;0;0],1,N,M);
diffInv = quatmult(q,qinv) - qId;

thr = 1e-10;
ok = all(abs(diffConj(:))<thr) && all(abs(diffDouble(:))<thr) && all(abs(diffInv(:))<thr);
