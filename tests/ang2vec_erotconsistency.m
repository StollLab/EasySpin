function ok = test()

% Check consistency of ang2vec() with with erot()

N = 500;
rng(88);
phi = rand(N,1)*2*pi;
theta = rand(N,1)*2*pi;

for k = 1:N
  n1 = ang2vec(phi(k),theta(k));
  R = erot(0,-theta(k),-phi(k));
  n2 = R(:,3);
  ok(k) = areequal(n1,n2,1e-10,'abs');
end
