function ok = test(opt)

%======================================================
% Two-pulse ESEEM
%======================================================
N = 200; dt = .01; %us
lw = 10; nlw = 1.57; cent = 0; noff = 101;
A = 7; B = 2; nuN = 14.9;
sys = [.5 .5];

[Sx,Sy,Sz] = sop(sys,'x1','y1','z1');
[Ix,Iz] = sop(sys,'x2','z2');

H = A*Sz*Iz + B*Sz*Ix + nuN*Iz;
Detector = Sx+1i*Sy;

td = zeros(N,1);
offset = nlw*lw*linspace(-1,1,noff);
weight = lshape(offset,cent,lw);
weight = weight/sum(weight);
  
P90 = expm(-1i*pi/2*Sy);
P180 = P90*P90;
prepDens = P90*(-Sz)*P90';

for iOff = 1:length(offset)
  td = td + weight(iOff)*evolve(prepDens, Detector, H+offset(iOff)*Sz, N,dt, [1 1], P180);
end

t = (0:N-1)*dt;
f = fdaxis(dt,N);
ttd = td - mean(td);
ttd = apowin('ham',N).*ttd;
fd = fftshift(abs(fft(ttd)));

ok = true;

if opt.Display
  subplot(2,2,1);
  plot(t,real(td),'r');
  axis tight
  xlabel('tau (us)');
  subplot(2,2,2);
  plot(t,imag(td),'b');
  axis tight
  xlabel('tau (us)');
  subplot(2,1,2);
  plot(f,real(fd));
  axis tight
  xlim([0 inf]);
  xlabel('frequency (MHz)');
end
