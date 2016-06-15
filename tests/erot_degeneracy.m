function [err,data] = test(opt,olddata)

% Test 3: alpha gamma degeneracy
%=========================================================
R0 = erot(pi/2,0,0);
R1 = erot(0,0,pi/2);
R2 = erot(pi/4,0,pi/4);

clear ok;
ok(1) = all(abs(R0(:)-R1(:))<1e-10);
ok(2) = all(abs(R0(:)-R2(:))<1e-10);

err = any(~ok);
data =[];
