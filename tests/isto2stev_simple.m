function ok = test()

% Test correctness of transformation matrix for all ranks

for k = 1:12
  n = 2*k+1;
  
  S = k/2; % minimum S needed for rank-k
  
  % Get transformation matrix
  C = isto2stev(k);
  
  % Get both operator sets
  T = cell(n,1);
  O = cell(n,1);
  for q = -k:k
    T{k+1-q} = isto(S,[k q]);
    O{k+1-q} = stev(S,[k,q]);
  end
  
  % Apply transformation to get O from T
  OO = cell(n,1);
  for o = 1:n
    OO{o} = 0;
    for t = 1:n
      OO{o} = OO{o} + C(o,t)*T{t};
    end
  end
  
  ok(k) = true;
  for o = 1:n
    ok(k) = ok(k) && areequal(OO{o},O{o},1e-10,'rel');
  end
  
end
