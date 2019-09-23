function [err,data] = test(opt,olddata)

%======================================================
% Simple axis permutations
%======================================================
thr = 1e-10;

R = [1 0 0; 0 0 1; 0 -1 0];
[n,phi] = rotmat2axi(R);
ok(1) = areequal(n,[1; 0; 0],thr,'abs') & areequal(phi,pi/2,thr,'abs');

R = [0 0 -1; 0 1 0; 1 0 0];
[n,phi] = rotmat2axi(R);
ok(2) = areequal(n,[0; 1; 0],thr,'abs') & areequal(phi,pi/2,thr,'abs');

R = [0 1 0; -1 0 0; 0 0 1];
[n,phi] = rotmat2axi(R);
ok(3) = areequal(n,[0; 0; 1],thr,'abs') & areequal(phi,pi/2,thr,'abs');

R = [0 1 0; 0 0 1; 1 0 0];
[n,phi] = rotmat2axi(R);
ok(4) = areequal(n,[1; 1; 1]/sqrt(3),thr,'abs') & areequal(phi,2*pi/3,thr,'abs');

err = any(~ok);
data = [];
