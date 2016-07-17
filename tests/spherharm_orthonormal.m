function [err,data] = test(opt,olddata)

% Test whether spherical and tesseral harmonics are orthonormal
%==========================================================================

thr = 1e-7;
fcn1 = @(theta,phi)...
  sin(theta).*conj(spherharm(3,2,theta,phi)).*spherharm(3,2,theta,phi);
v = dblquad(fcn1,0,pi,0,2*pi,thr,@quadl);
err(1) = abs(v-1)>10*thr;

fcn2 = @(theta,phi)...
  sin(theta).*conj(spherharm(2,2,theta,phi)).*spherharm(2,1,theta,phi);
v = dblquad(fcn2,0,pi,0,2*pi,thr,@quadl);
err(2) = (abs(v)>10*thr);

fcn3c = @(theta,phi)...
  sin(theta).*spherharm(3,1,theta,phi,'c').*spherharm(3,1,theta,phi,'c');
v = dblquad(fcn3c,0,pi,0,2*pi,thr,@quadl);
err(3) = (abs(v-1)>10*thr);

fcn3s = @(theta,phi)...
  sin(theta).*spherharm(3,1,theta,phi,'s').*spherharm(3,1,theta,phi,'s');
v = dblquad(fcn3s,0,pi,0,2*pi,thr,@quadl);
err(4) = (abs(v-1)>10*thr);

err = any(err);

data = [];
