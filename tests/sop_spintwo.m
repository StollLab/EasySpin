function ok = test()

% Check matrices of all basic spin operators for S=2

S = 2;

ops = 'exyz+-';
for k = 1:numel(ops)
  Op{k} = sop(S,ops(k));
end

a = sqrt(3/2);
b = sqrt(6);

myOp{1} = eye(5);      % identity
myOp{2} = [0 1 0 0 0; 1 0 a 0 0; 0 a 0 a 0; 0 0 a 0 1; 0 0 0 1 0];  % Sx
myOp{3} = [0 1 0 0 0; -1 0 a 0 0; 0 -a 0 a 0; 0 0 -a 0 1; 0 0 0 -1 0]/1i;  % Sy
myOp{4} = diag(2:-1:-2);  % Sz
myOp{5} = diag([2 b b 2],+1);  % S+
myOp{6} = diag([2 b b 2],-1);  % S-

for k = 1:numel(Op)
  ok(k) = areequal(Op{k},myOp{k},1e-12,'abs');
end
