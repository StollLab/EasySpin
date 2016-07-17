function [err,data] = test(opt,olddata)

% Test 2: S=1
%======================================================
Q{1} = [1 1 1];
Q{2} = [1 2 3 2 1];
Q{3} = [1 3 6 7 6 3 1];
Q{4} = [1 4 10 16 19 16 10 4 1];

clear ok
for n = 1:numel(Q)
  q = equivsplit(1,n);
  ok(n) = all(q==Q{n});
end

err = any(~ok);
data = [];
