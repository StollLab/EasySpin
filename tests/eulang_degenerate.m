function ok = test()

% For beta = 0 and beta = pi, alpha and gamma cannot be separated. eulang
% sets gamma to zero and puts the entire z rotation into alpha. This is
% alpha+gamma for beta = 0 and alpha-gamma for beta = pi.

alpha = 0.5;
gamma = 0.3;

% beta = 0
ang = eulang(erot(alpha,0,gamma));
ok(1) = areequal(ang,[alpha+gamma 0 0],1e-12,'abs');

% beta = pi
ang = eulang(erot(alpha,pi,gamma));
ok(2) = areequal(ang,[alpha-gamma pi 0],1e-12,'abs');

% beta = pi, with alpha-gamma negative: alpha must be shifted into [0,2*pi)
ang = eulang(erot(gamma,pi,alpha));
ok(3) = areequal(ang,[gamma-alpha+2*pi pi 0],1e-12,'abs');
