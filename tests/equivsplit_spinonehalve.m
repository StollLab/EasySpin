function [err,data] = test(opt,olddata)

% Test 1: S=1/2
%======================================================
Q{1} = [1 1];
Q{2} = [1 2 1];
Q{3} = [1 3 3 1];
Q{4} = [1 4 6 4 1];
Q{5} = [1 5 10 10 5 1];

clear ok
for n = 1:numel(Q)
  q = equivsplit(1/2,n);
  ok(n) = all(q==Q{n});
end

err = any(~ok);

data = [];
