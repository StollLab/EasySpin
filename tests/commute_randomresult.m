function [err,data] = test(opt,olddata)

N = 10;
A = complex(rand(N),rand(N));
B = complex(rand(N),rand(N));
r1 = commute(A,B);
r2 = A*B - B*A;
err = max(abs(r1(:)-r2(:))) > 1e-10;

data = [];
