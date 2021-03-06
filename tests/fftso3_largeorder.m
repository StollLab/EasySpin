function ok = test()

Lmax = 30;
LMK = [20 11 -5];
A = 1+2i;

f = @(a,b,c) A*wignerd(LMK,a,b,c);

[LMK_,A_] = fftso3(f,Lmax);

ok = size(LMK_,1)==1 && all(LMK_==LMK) && areequal(A,A_,1e-10,'abs');

