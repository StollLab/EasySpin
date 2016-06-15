function ok = areequal(A,B,delta)

if (nargin<3), delta = 1e-10; end

ok = all(abs(A(:)-B(:))<delta);

return
