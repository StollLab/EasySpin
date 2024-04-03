function ok = test(opt)

%=================================================================
% Test the use of multiple detection operators with evolve
%=================================================================
N = 100; dt = .005;  % µs
lw = 10; nlw = 1.57; cent = 0; noff = 101;
A = 7; B = 5; nuN = 14.9; % MHz
sys = [1/2 1/2];

[Sx,Sy,Sz] = sop(sys,'x1','y1','z1');
[Ix,Iz] = sop(sys,'x2','z2');

H = A*Sz*Iz + B*Sz*Ix + nuN*Iz;
Detectors = {Sx,Sy,Sz};

sigx = zeros(N,1);
sigy = zeros(N,1);
sigz = zeros(N,1);

offset = nlw*lw*linspace(-1,1,noff);
weight = lshape(offset,cent,lw);
weight = weight/sum(weight);
  
P90 = expm(-1i*pi/2*0.8*Sy);
P180 = P90*P90;
prepDens = P90*(-Sz)*P90';

for iOff = 1:length(offset)
  a = evolve(prepDens, Detectors, H+offset(iOff)*Sz, N,dt, [1 1], P180);
  sigx = sigx + weight(iOff)*a{1};
  sigy = sigy + weight(iOff)*a{2};
  sigz = sigz + weight(iOff)*a{3};
end
sigx = real(sigx);
sigy = real(sigy);
sigz = real(sigz);

t = (0:N-1)*dt;

ok = true;

if opt.Display
  subplot(1,1,1);
  plot(t,sigx,'r',t,sigy,'g',t,sigz,'b');
  axis tight
  xlabel('tau (µs)');
end
