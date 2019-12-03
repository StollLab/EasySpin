function [err,data] = test(opt,olddata)

% Check matrices of all spin operators for S=1/2

ops = 'exyz+-ab';
for k = 1:numel(ops)
  Op{k} = sop(1/2,ops(k));
end

myOp{1} = [1 0; 0 1];
myOp{2} = [0 1;1 0]/2;
myOp{3} = [0 -1;1 0]*1i/2;
myOp{4} = [1 0;0 -1]/2;
myOp{5} = [0 1;0 0];
myOp{6} = [0 0; 1 0];
myOp{7} = [1 0; 0 0];
myOp{8} = [0 0; 0 1];

for k = 1:numel(Op)
  ok(k) = areequal(Op{k},myOp{k},1e-12,'abs');
end

err = any(~ok);

data = [];

