function ok = test()

% Check matrices of all spin operators for S=1/2

ops = 'exyz+-ab';
for k = 1:numel(ops)
  Op{k} = sop(1/2,ops(k));
end

myOp{1} = [1 0; 0 1];      % identity
myOp{2} = [0 1;1 0]/2;     % Sx
myOp{3} = [0 -1;1 0]*1i/2; % Sy
myOp{4} = [1 0;0 -1]/2;    % Sz
myOp{5} = [0 1;0 0];       % S+
myOp{6} = [0 0; 1 0];      % S-
myOp{7} = [1 0; 0 0];      % Sa
myOp{8} = [0 0; 0 1];      % Sb

for k = 1:numel(Op)
  ok(k) = areequal(Op{k},myOp{k},1e-12,'abs');
end
