function ok = test()

% Check matrices of all basic spin operators for S=3/2

ops = 'exyz+-';
for k = 1:numel(ops)
  Op{k} = sop(3/2,ops(k));
end

a = sqrt(3/4);
b = sqrt(3);

myOp{1} = [1 0 0 0;  0 1 0 0; 0  0  1 0; 0 0  0  1];      % identity
myOp{2} = [0 a 0 0;  a 0 1 0; 0  1  0 a; 0 0  a  0];      % Sx
myOp{3} = [0 a 0 0; -a 0 1 0; 0 -1  0 a; 0 0 -a  0]/1i;   % Sy
myOp{4} = [3 0 0 0;  0 1 0 0; 0  0 -1 0; 0 0  0 -3]/2;    % Sz
myOp{5} = [0 b 0 0;  0 0 2 0; 0  0  0 b; 0 0  0  0];      % S+
myOp{6} = [0 0 0 0;  b 0 0 0; 0  2  0 0; 0 0  b  0];      % S-

for k = 1:numel(Op)
  ok(k) = areequal(Op{k},myOp{k},1e-12,'abs');
end
