function ok = test()

% Make sure sop() tolerates zero spin

S = [1/2 0];
Op = sop(S,'x1');
ok(1) = length(Op)==prod(2*S+1);

S = [0 1/2];
Op = sop(S,'x2');
ok(2) = length(Op)==prod(2*S+1);

S = [1/2 0 1/2];
Op = sop(S,'x1y3');
ok(3) = length(Op)==prod(2*S+1);
