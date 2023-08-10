function ok = test()

% Cubic high-spin systems

sys.S = 7/2;
sys.g = 2*[1 1 1];
B2 = 0;
B4 = 0;
B6 = 10;

sys.B2 = B2;
sys.B4 = [5*B4 0 0 0 B4 0 0 0 0];
sys.B6 = [0 0 -21*B6 0 0 0 B6 0 0 0 0 0 0];
G = hamsymm(sys);
ok(1) = strcmp(G,'Oh');

B2 = 0;
B4 = 10;
B6 = 0;
sys.B2 = B2;
sys.B4 = [5*B4 0 0 0 B4 0 0 0 0];
sys.B6 = [0 0 -21*B6 0 0 0 B6 0 0 0 0 0 0];
G = hamsymm(sys);
ok(2) = strcmp(G,'Oh');

B2 = 0;
B4 = 10;
B6 = 12.87;
sys.B2 = B2;
sys.B4 = [5*B4 0 0 0 B4 0 0 0 0];
sys.B6 = [0 0 -21*B6 0 0 0 B6 0 0 0 0 0 0];

G = hamsymm(sys);
ok(3) = strcmp(G,'Oh');
