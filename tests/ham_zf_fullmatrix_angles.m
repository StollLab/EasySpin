function ok = test()

% full D matrices with Euler angles

Sys.S = 3/2;

D = rand(3);
DFrame = rand(1,3)*2*pi;
Sys.D = D;
Sys.DFrame = DFrame;

H1 = ham_zf(Sys);

R = erot(DFrame);
Dr = R.'*D*R;
Sys.D = Dr;
Sys.DFrame = [0 0 0];
H2 = ham_zf(Sys);

ok = areequal(H1,H2,1e-7,'rel');
