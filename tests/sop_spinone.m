function ok = test()

% Check all matrices for S=1

ops = 'exyz+-';
for k = 1:numel(ops)
  Op{k} = sop(1,ops(k));
end

myOp{1} = [1 0 0; 0 1 0; 0 0 1];
myOp{2} = [0 1 0;1 0 1; 0 1 0]/sqrt(2);
myOp{3} = [0 -1 0;1 0 -1; 0 1 0]*1i/sqrt(2);
myOp{4} = [1 0 0; 0 0 0;0 0 -1];
myOp{5} = [0 1 0;0 0 1; 0 0 0]*sqrt(2);
myOp{6} = [0 0 0; 1 0 0; 0 1 0]*sqrt(2);

for k = 1:numel(Op)
  ok(k) = areequal(Op{k},myOp{k},1e-12,'abs');
end
