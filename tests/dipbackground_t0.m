function ok = test()

% At t=0, there is no decay yet, so Vinter must equal 1 regardless of
% concentration and modulation depth.

conc = [10 100 500];   % µM
lambda = [0.1 0.4 0.9];

ok = true(numel(conc),numel(lambda));
for i = 1:numel(conc)
  for j = 1:numel(lambda)
    Vinter0 = dipbackground(0,conc(i),lambda(j));
    ok(i,j) = areequal(Vinter0,1,1e-12,'abs');
  end
end

end
