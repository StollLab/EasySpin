function ok = test()

% Explicit check of transformation matrix for S1=S2=1

U2C = cgmatrix(1,1);

A = 1/sqrt(2);
B = 1/sqrt(6);
C = sqrt(2/3);
D = 1/sqrt(3);

U2c_ref = [
  1  0  0  0  0  0  0  0  0;
  0  A  0  A  0  0  0  0  0;
  0  0  B  0  C  0  B  0  0;
  0  0  0  0  0  A  0  A  0;
  0  0  0  0  0  0  0  0  1;
  0  A  0 -A  0  0  0  0  0;
  0  0  A  0  0  0 -A  0  0;
  0  0  0  0  0  A  0 -A  0;
  0  0  D  0 -D  0  D  0  0];

ok = areequal(U2C,U2c_ref,1e-10,'abs');
