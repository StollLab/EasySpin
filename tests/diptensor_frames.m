function ok = test()

rng(34);

r = (2*rand(1,3)-1)*5;  % distance between two spins, nm

% g tensor principal values and Euler angles for both spins
gpv1 = [2 2.1 2.2];
gFrame1 = rand(1,3)*pi;
gpv2 = [2 2.1 2.2];
gFrame2 = rand(1,3)*pi;

% Calculate full g tensors
R1 = erot(gFrame1).';
g1 = R1*diag(gpv1)*R1.';
R2 = erot(gFrame2).';
g2 = R2*diag(gpv2)*R2.';

T_full = diptensor(g1,g2,r);
T_pvframes = diptensor({gpv1,gFrame1},{gpv2,gFrame2},r);

ok = areequal(T_full,T_pvframes,1e-10,'abs');
