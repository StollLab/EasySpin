function ok = test()

rho = 2*pi/5;  % rotation angle, radians

angles0 = [pi/5, pi/6, pi/7];

axletter = {'x','y','z','xy','xz','yz','xyz'};
axvector = {[1 0 0],[0 1 0],[0 0 1],[1 1 0],[1 0 1],[0 1 1],[1 1 1]};

for k = numel(axletter):-1:1
  angles1 = rotateframe(angles0,axletter{k},rho);
  angles2 = rotateframe(angles0,axvector{k},rho);
  ok(k) = areequal(angles1,angles2,1e-10,'abs');
end
