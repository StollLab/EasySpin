function ok = test()

% eulang must return alpha and gamma in [0,2*pi) and beta in [0,pi], also
% for rotation matrices generated from angles outside these ranges.

rng(2342);

% Example from the documentation: gamma = -143 deg is returned as 217 deg
ang = eulang(erot([34 72 -143]*pi/180))*180/pi;
ok(1) = areequal(ang,[34 72 217],1e-10,'abs');

% Random angles, including negative ones and ones larger than 2*pi
ang0 = (rand(20,3)-0.5)*6*pi;
for k = size(ang0,1):-1:1
  R = erot(ang0(k,:));
  ang = eulang(R);
  inRange = ang(1)>=0 && ang(1)<2*pi && ...
            ang(2)>=0 && ang(2)<=pi && ...
            ang(3)>=0 && ang(3)<2*pi;
  ok(k+1) = inRange && areequal(erot(ang),R,1e-12,'abs');
end
