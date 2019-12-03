function [err,data] = test(opt,olddata)

% Test S_alpha and S_beta for spin 1/2
%-------------------------------------------------------

Op{1} = sop(1/2,'a');
Op{2} = sop(1/2,'b');

myOp{1} = [1 0; 0 0];
myOp{2} = [0 0; 0 1];

for k = 1:numel(Op)
  ok(k) = areequal(Op{k},myOp{k},1e-10,'abs');
end

err = any(~ok);

data = [];

