function ok = test()

% S=1/2 with multiple nuclear spins-2

Q{1} = [1 1 1 1 1];
Q{2} = [1 2 3 4 5 4 3 2 1];
Q{3} = [1 3 6 10 15 18 19 18 15 10 6 3 1];
Q{4} = [1 4 10 20 35 52 68 80 85 80 68 52 35 20 10 4 1];
Q{5} = [1 5 15 35 70 121 185 255 320 365 381 365 320 255 185 121 70 35 15 5 1];

for n = 1:numel(Q)
  q = equivsplit(2,n);
  ok(n) = numel(q)==numel(Q{n}) && all(q==Q{n});
end
