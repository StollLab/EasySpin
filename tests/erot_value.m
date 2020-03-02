function ok = test()

alpha = 1.3*pi;
beta = 0.4*pi;
gamma = -0.8*pi;

Rg = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
Rb = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
Ra = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
R0 = Rg*Rb*Ra;

R1 = erot(alpha,beta,gamma);
R2 = erot([alpha,beta,gamma]);

ok = areequal(R1,R0,1e-10,'abs') && areequal(R2,R0,1e-10,'abs');
