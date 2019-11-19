function err = test()

% Comparison of two output modes: full matrix vs. rows
%======================================================
a = [1.2 0.654 0.912334]*pi;

R1 = erot(a);
[x,y,z] = erot(a,'rows');
R2 = [x y z].';

err = ~areequal(R1,R2);
