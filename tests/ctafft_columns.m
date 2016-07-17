function [err,data] = test(opt,olddata)

%======================================================
% Test 4: Check whether ctafft operates along columns
%======================================================
TD = rand(512,2);
TD(:,3) = 0;
FD = ctafft(TD,1:5);
err = any(FD(:,3)~=0);
data = [];
