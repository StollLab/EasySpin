function ok = test()

% Multiple spins-3/2

test(1).n = 2;
test(1).F = [3 2 1 0];
test(1).N = [1 1 1 1];

test(2).n = 3;
test(2).F = [4.5 3.5 2.5 1.5 0.5];
test(2).N = [1 2 3 4 2];

test(3).n = 4;
test(3).F = [6 5 4 3 2 1 0];
test(3).N = [1 3 6 10 11 9 4];

test(4).n = 5;
test(4).F = [7.5 6.5 5.5 4.5 3.5 2.5 1.5 0.5];
test(4).N = [1 4 10 20 30 36 34 20];

S = 3/2;
for q=1:numel(test)
  [F,N] = equivcouple(S,test(q).n);
  ok(q) = all(F==test(q).F) & all(N==test(q).N);
  dim = (2*S+1)^test(q).n;
  ok(q) = ok(q) && (sum((2*F+1).*N)==dim);
end
