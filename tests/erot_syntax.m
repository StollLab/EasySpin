function [err,data] = test(opt,olddata)

% Test 1: syntax
%======================================================
a = rand(1,3);
R1 = erot(a);
R2 = erot(a(1),a(2),a(3));
err = any(R1(:)~=R2(:));

data = [];
