function [err,data] = test(opt,olddata)

% Test 2: all matrices for S=1

ops = 'exyz+-';
for k=1:numel(ops);
  Op{k} = sop(1,ops(k));
end
myOp{1} = [1 0 0; 0 1 0; 0 0 1];
myOp{2} = [0 1 0;1 0 1; 0 1 0]/sqrt(2);
myOp{3} = [0 -1 0;1 0 -1; 0 1 0]*1i/sqrt(2);
myOp{4} = [1 0 0; 0 0 0;0 0 -1];
myOp{5} = [0 1 0;0 0 1; 0 0 0]*sqrt(2);
myOp{6} = [0 0 0; 1 0 0; 0 1 0]*sqrt(2);

for k=1:numel(Op)
  ok(k) = all(abs(myOp{k}(:)-Op{k}(:))<1e-10);
end

if any(~ok)
  err = 1;
else
  err = 0;
end

data = [];

