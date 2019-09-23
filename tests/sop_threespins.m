function [err,data] = test(opt,olddata)

% Test syntax for three spins
%======================================================
kron3 = @(s,A,B,C)kron(kron(sop(s(1),A),sop(s(2),B)),sop(s(3),C));

s = [3/2 1/2 1];
Op{1} = sop(s,'x1y3');
myOp{1} = kron3(s,'x','e','y');

Op{2} = sop(s,'+2');
myOp{2} = kron3(s,'e','+','e');

Op{3} = sop(s,'z1z2z3');
myOp{3} = kron3(s,'z','z','z');

Op{4} = sop(s,'x3y1');
myOp{4} = kron3(s,'y','e','x');

for k = 1:numel(Op)
  ok(k) = areequal(Op{k},myOp{k},1e-10,'abs');
end

err = any(~ok);

data = [];
