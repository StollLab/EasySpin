function ok = test()

% Check input and output syntax, output size

S1 = 2;
S2 = 3/2;
N = hsdim([S1 S2]);
C = cgmatrix(S1,S2);

ok = numel(C)==N^2;

