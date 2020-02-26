function ok = test()

% Check all single-transition matrices for a single spin with S>1/2
%-------------------------------------------------------------------------------
S = 3/2;
components = 'xyz+-';
nLevels = 2*S+1;

idx = 1;
threshold = 1e-14;

% Compare diagonal single-element operators: 'e(L)' and 'e(L|L)'
for k = 1:numel(components)
  for L = 1:nLevels
    ops1 = sprintf('e(%d)',L);
    ops2 = sprintf('e(%d|%d)',L,L);
    sopOp1 = sop(S,ops1);
    sopOp2 = sop(S,ops2);
    correctOp = zeros(2*S+1);
    correctOp(L,L) = 1;
    ok(idx) = areequal(sopOp1,correctOp,threshold,'abs');
    ok(idx+1) = areequal(sopOp2,correctOp,threshold,'abs');
    idx = idx+2;
  end
end

% Compare transition-selective x,y,z,+,- operators: 'x(L1|L2)' etc.
for k = 1:numel(components)
  for L1 = 1:nLevels
    for L2 = L1+1:nLevels
      ops = sprintf('%s(%d|%d)',components(k),L1,L2);
      sopOp = sop(S,ops);
      correctOp = zeros(nLevels);
      switch components(k)
        case 'x', correctOp(L1,L2) = 1/2; correctOp(L2,L1) = 1/2;
        case 'y', correctOp(L1,L2) = 1/2i; correctOp(L2,L1) = -1/2i;
        case 'z', correctOp(L1,L1) = 1/2; correctOp(L2,L2) = -1/2;
        case '+', correctOp(L1,L2) = 1;
        case '-', correctOp(L2,L1) = 1;
      end
      ok(idx) = areequal(sopOp,correctOp,threshold,'abs');
      idx = idx + 1;
    end
  end
end
