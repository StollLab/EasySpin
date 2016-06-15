function [err,data] = test(opt,olddata)

% All matrices for S=1/2

ops = 'exyz+-';
for k=1:numel(ops)
  Op{k} = sop(1/2,ops(k));
end
myOp{1} = [1 0; 0 1];
myOp{2} = [0 1;1 0]/2;
myOp{3} = [0 -1;1 0]*1i/2;
myOp{4} = [1 0;0 -1]/2;
myOp{5} = [0 1;0 0];
myOp{6} = [0 0; 1 0];

for k=1:numel(Op)
  ok(k) = all(myOp{k}(:)==Op{k}(:));
end

err = any(~ok);

data = [];

