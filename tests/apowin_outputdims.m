function [err,data] = test(opt,olddata)

%======================================================
% Test 3: test dimension of output
%======================================================

N = 2345;
w = apowin('ham',N);
err = (size(w,1)~=N) | (size(w,2)~=1);

data = [];
