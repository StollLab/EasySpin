function ok = test()

% S=1/2 with multiple nuclear spins-3/2

Q{1} = [1 1 1 1];
Q{2} = [1 2 3 4 3 2 1];
Q{3} = [1 3 6 10 12 12 10 6 3 1];
Q{4} = [1 4 10 20 31 40 44 40 31 20 10 4 1];
Q{5} = [1 5 15 35 65 101 135 155 155 135 101 65 35 15 5 1];

clear ok
for n = 1:numel(Q)
  q = equivsplit(3/2,n);
  ok(n) = numel(q)==numel(Q{n}) && all(q==Q{n});
end
