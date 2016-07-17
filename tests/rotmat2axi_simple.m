function [err,data] = test(opt,olddata)

%======================================================
% Simple axis permutations
%======================================================
R = [1 0 0; 0 0 1; 0 -1 0];
[n,phi] = rotmat2axi(R);
ok(1) = areequal(n,[1 0 0]) & areequal(phi,pi/2);

R = [0 0 -1; 0 1 0; 1 0 0];
[n,phi] = rotmat2axi(R);
ok(2) = areequal(n,[0 1 0]) & areequal(phi,pi/2);

R = [0 1 0; -1 0 0; 0 0 1];
[n,phi] = rotmat2axi(R);
ok(3) = areequal(n,[0 0 1]) & areequal(phi,pi/2);

R = [0 1 0; 0 0 1; 1 0 0];
[n,phi] = rotmat2axi(R);
ok(4) = areequal(n,[1 1 1]/sqrt(3)) & areequal(phi,2*pi/3);

err = any(~ok);
data = [];
