function ok = test()

% Test whether spherical and tesseral harmonics are orthonormal

thr = 1e-7;
fcn1 = @(theta,phi)...
  sin(theta).*conj(spherharm(3,2,theta,phi)).*spherharm(3,2,theta,phi);
v = integral2(fcn1,0,pi,0,2*pi,'AbsTol',thr);
ok(1) = abs(v-1)<10*thr;

fcn2 = @(theta,phi)...
  sin(theta).*conj(spherharm(2,2,theta,phi)).*spherharm(2,1,theta,phi);
v = integral2(fcn2,0,pi,0,2*pi,'AbsTol',thr);
ok(2) = (abs(v)<10*thr);

fcn3c = @(theta,phi)...
  sin(theta).*spherharm(3,1,theta,phi,'c').*spherharm(3,1,theta,phi,'c');
v = integral2(fcn3c,0,pi,0,2*pi,'AbsTol',thr);
ok(3) = (abs(v-1)<10*thr);

fcn3s = @(theta,phi)...
  sin(theta).*spherharm(3,1,theta,phi,'s').*spherharm(3,1,theta,phi,'s');
v = integral2(fcn3s,0,pi,0,2*pi,'AbsTol',thr);
ok(4) = (abs(v-1)<10*thr);
