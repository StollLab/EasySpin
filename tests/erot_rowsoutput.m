function [err,data] = test(opt,olddata)

% Comparison of two output modes: full matrix vs. rows
%======================================================
a = rand(1,3);
R1 = erot(a);
[x,y,z] = erot(a,'rows');
R2 = [x y z].';
err = any(R1(:)~=R2(:));

data = [];
