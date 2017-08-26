function [err,data] = test(opt,olddata)

% Check all matrices for S=1

ops = {'x(1|2)' 'p(1|2)1' 'y(2|3)' 'z(2|3)' 'e(3|3)' 'e(2)'};

for k = 1:numel(ops)
  Op{k} = sop(1,ops{k});
end
myOp{1} = [0 0.5 0; 0.5 0 0; 0 0 0];
myOp{2} = [0 1 0;0 0 0; 0 0 0];
myOp{3} = [0 0 0;0 0 -0.5; 0 0.5 0]*1i;
myOp{4} = [0 0 0; 0 1/2 0;0 0 -1/2];
myOp{5} = [0 0 0;0 0 0; 0 0 1];
myOp{6} = [0 0 0; 0 1 0; 0 0 0];

for k = 1:numel(Op)
  ok(k) = areequal(Op{k},myOp{k},1e-12);
end

err = any(~ok);

data = [];

