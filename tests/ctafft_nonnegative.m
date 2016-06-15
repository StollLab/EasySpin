function [err,data] = test(opt,olddata)

%======================================================
% Test 2: Output should be nonnegative
%======================================================
TD = rand(1,1024);
FD = ctafft(TD,1:5,2048);
err = ~isreal(FD) | any(FD(:)<0);

data = [];
