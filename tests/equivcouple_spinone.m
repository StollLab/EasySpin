function ok = test()

% Multiple spins-1

test(1).n = 2;
test(1).F = [2 1 0];
test(1).N = [1 1 1];

test(2).n = 3;
test(2).F = [3 2 1 0];
test(2).N = [1 2 3 1];

test(3).n = 4;
test(3).F = [4 3 2 1 0];
test(3).N = [1 3 6 6 3];

test(4).n = 5;
test(4).F = [5 4 3 2 1 0];
test(4).N = [1 4 10 15 15 6];

test(5).n = 6;
test(5).F = [6 5 4 3 2 1 0];
test(5).N = [1 5 15 29 40 36 15];

S = 1;
for q=1:numel(test)
  [F,N] = equivcouple(S,test(q).n);
  ok(q) = all(F==test(q).F) & all(N==test(q).N);
  dim = (2*S+1)^test(q).n;
  ok(q) = ok(q) && (sum((2*F+1).*N)==dim);
end
