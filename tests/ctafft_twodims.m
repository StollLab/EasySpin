function [err,data] = test(opt,olddata)

%======================================================
% Test 3: Check output size for 2D case
%======================================================
TD = rand(512,20);
FD = ctafft(TD,1:5);
err = any(size(FD)~=[512 20]);
data = [];
