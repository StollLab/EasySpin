function [err,data] = test(opt,olddata)

% Test some matrices for S=I=1/2
%======================================================
s = [1/2 1/2];
Op{1} = sop(s,'xx');
myOp{1} = [0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0]/4;

Op{2} = sop(s,'z+');
myOp{2} = [0 1 0 0; 0 0 0 0; 0 0 0 -1; 0 0 0 0]/2;

Op{3} = sop(s,'ey');
myOp{3} = [0 -1 0 0; 1 0 0 0; 0 0 0 -1; 0 0 1 0]*1i/2;

Op{4} = sop(s,'--');
myOp{4} = [0 0 0 0; 0 0 0 0; 0 0 0 0; 1 0 0 0];

for k = 1:numel(Op)
  ok(k) = areequal(Op{k},myOp{k},1e-12,'abs');
end

err = any(~ok);

data = [];
