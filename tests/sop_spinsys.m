function [err,data] = test(opt,olddata)

% Assert that sop can handle spin systems as first input
%-------------------------------------------------------

Sys{1}.S = 1/2;
Sys{2}.S = [3/2 1]; Sys{2}.ee = 1;
Sys{3}.S = 1;
Sys{3}.Nucs = '63Cu'; Sys{3}.A = 1;

for k = 1:numel(Sys)
  Spins = spinvec(Sys{k});
  op = repmat('x',1,numel(Spins));
  Op1 = sop(Sys{k},op);
  Op2 = sop(Spins,op);
  ok(k) = areequal(Op1,Op2,1e-12,'abs');
end

err = ~ok;

data = [];
