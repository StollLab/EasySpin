function ok = test()

% Check some transition-selective operators for a two-spin system

Spins = [1/2 1];
ops = {'x(1|2)1' '+(1|2),e(1)' 'x,z' 'z2,x1' ...
       'y(1|2)1,+(1|2)2' 'e(1)1,e(2)2' 'x1,x(1|3)2' ...
       'b,-(1|3)'};

for k = 1:numel(ops)
  Op{k} = sop(Spins,ops{k});
end

myOp{1} = kron([0 0.5; 0.5 0],[1 0 0; 0 1 0; 0 0 1]);
myOp{2} = kron([0 1; 0 0],[1 0 0; 0 0 0; 0 0 0]);
myOp{3} = kron([0 0.5; 0.5 0],[1 0 0; 0 0 0; 0 0 -1]);
myOp{4} = myOp{3};
myOp{5} = kron([0 -0.5i; 0.5i 0],[0 1 0; 0 0 0; 0 0 0]);
myOp{6} = kron([1 0; 0 0],[0 0 0; 0 1 0; 0 0 0]);
myOp{7} = kron([0 0.5; 0.5 0],[0 0 0.5; 0 0 0; 0.5 0 0]);
myOp{8} = kron([0 0; 0 1],[0 0 0; 0 0 0; 1 0 0]);

for k = 1:numel(Op)
  ok(k) = areequal(Op{k},myOp{k},1e-12,'abs');
end
