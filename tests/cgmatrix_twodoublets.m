function [err,data] = test(opt,olddata)

% Explicit check of transformation matrix for S1=S2=1/2
%======================================================

U2C = cgmatrix(1/2,1/2);

c = 1/sqrt(2);
U2c_ref = [1 0 0 0; 0 c c 0; 0 0 0 1; 0 c -c 0];

err = ~areequal(U2C,U2c_ref,1e-10,'abs');
data = [];
