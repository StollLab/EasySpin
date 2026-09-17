function ok = test()

% K = dipkernel(t,r) must be a numel(t) x numel(r) matrix, regardless
% of whether t and r are given as row or column vectors.

t = linspace(-1,2,15);  % µs
r = linspace(2,5,10);   % nm

K = dipkernel(t,r);
ok(1) = isequal(size(K),[numel(t) numel(r)]);

Kc = dipkernel(t(:),r(:));
ok(2) = isequal(size(Kc),[numel(t) numel(r)]);
ok(3) = areequal(K,Kc,1e-12,'abs');

end
