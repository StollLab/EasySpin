function ok = test()

% full hyperfine matrices with Euler angles

Sys.S = 1/2;
Sys.Nucs = '1H,14N';
A1 = rand(3);
A2 = rand(3);
AFrame1 = rand(1,3)*pi;
AFrame2 = rand(1,3)*pi;
Sys.A = [A1; A2];
Sys.AFrame = [AFrame1; AFrame2];
H1 = hfine(Sys);

R1 = erot(AFrame1);
R2 = erot(AFrame2);
A1r = R1.'*A1*R1;
A2r = R2.'*A2*R2;
Sys.A = [A1r; A2r];
Sys.AFrame = [0 0 0; 0 0 0];
H2 = hfine(Sys);

deviation = max(abs(H1(:)-H2(:)))/max(abs(H1(:)));

ok = deviation<1e-7;
