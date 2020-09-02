function ok = test()

% Test whether spherical harmonics are normalized to 1

LM = [0 0; 1 -1; 1 0; 1 1; 2 -2; 2 -1; 2 0; 2 1; 2 2; 3 -3; 3 -2; 3 -1; 3 0; 3 1; 3 2; 3 3];

thr = 1e-6;
for f = 1:size(LM,1)
  L = LM(f,1);
  M = LM(f,2); 
  fcn = @(theta,phi) sin(theta).*conj(spherharm(L,M,theta,phi,'r')).*spherharm(L,M,theta,phi,'r');
  v = real(integral2(fcn,0,pi,0,2*pi,'RelTol',thr));
  ok(f) = areequal(v,1,thr,'rel');
end
