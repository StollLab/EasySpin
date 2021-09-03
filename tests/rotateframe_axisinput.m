function ok = test()

rho = 2*pi/5;  % rotation angle, radians

Ori = [pi/5, pi/6, pi/7];

angles1 = rotateframe(Ori,'x',rho);
angles2 = rotateframe(Ori,[1,0,0],rho);

ok = areequal(angles1,angles2,1e-10,'abs');
