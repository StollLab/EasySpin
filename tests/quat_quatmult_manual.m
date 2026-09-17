function ok = test()

% Check quatmult() against an independent, element-by-element
% implementation of the Hamilton product (using the 4x4 matrix
% representation of quaternion multiplication), for a multi-dimensional
% array of quaternions.

N = 5;
M = 10;

q = rand(4,N,M);
q = q./sqrt(sum(q.*q,1));

r = rand(4,N,M);
r = r./sqrt(sum(r.*r,1));

t_manual = zeros(4,N,M);
for i = 1:N
  for j = 1:M
    r0 = r(1,i,j); r1 = r(2,i,j); r2 = r(3,i,j); r3 = r(4,i,j);
    Mr = [r0 -r1 -r2 -r3;
          r1  r0  r3 -r2;
          r2 -r3  r0  r1;
          r3  r2 -r1  r0];
    t_manual(:,i,j) = Mr*q(:,i,j);
  end
end

t = quatmult(q,r);

diff = t - t_manual;

thr = 1e-10;
ok = all(abs(diff(:))<thr);
