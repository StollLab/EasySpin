function [err,data] = test(opt,olddata)

% S_alpha and S_beta for spin 1/2

Op{1} = sop(1/2,'a');
Op{2} = sop(1/2,'b');

myOp{1} = [1 0; 0 0];
myOp{2} = [0 0; 0 1];

for k=1:numel(Op)
  ok(k) = all(abs(myOp{k}(:)-Op{k}(:))<1e-10);
end

if any(~ok)
  err = 1;
else
  err = 0;
end

data = [];

