function ok = test()

% rotaxi2mat should the inverse of rotmat2axi

rng(44);

Angles = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; ...
  1 1 1; 0.5 0 0; 0 0.5 0; 0 0 0.5; 0.5 0.5 0; 0.5 0 0.5; ...
  0 0.5 0.5; 0.5 0.5 0.5; randn(10,3)]*pi;

for k = 1:size(Angles,1)
  R1 = erot(Angles(k,:));
  [n,phi] = rotmat2axi(R1);
  R2 = rotaxi2mat(n,phi);
  ok(k) = areequal(R1,R2,1e-8,'abs');
end
