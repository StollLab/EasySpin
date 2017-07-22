% excitation profile of a mw pulse
%==========================================================================
% In this example, we compute the excitation profiles of a pi/2 and a pi
% pulse of a given length and compare them to the sinc function (Fourier
% transform of the rectangular pulse shape), which we would expect
% if spin dynamics were linear.

clear, clf

% parameters
%-------------------------------------------------------------
tp = 0.064; % pulse length, µs
maxoffset = 3/tp; % offset range limit, MHz

% computation
%-------------------------------------------------------------
n = 4e3; % number of points
M0 = [0;0;1]; % equilibrium magnetization along z axis
offsetFreq = maxoffset*linspace(-1,1,n);
offPhase = 2*pi*tp*offsetFreq;

% loop over all offsets and compute z magnetization
flipAngle = pi;
Mz180 = zeros(1,n);
for k = 1:length(offsetFreq)
  M = expm([0 offPhase(k) 0; -offPhase(k) 0 -flipAngle; 0 flipAngle 0])*M0;
  Mz180(k) = M(3);
end

flipAngle = pi/2;
Mz90 = zeros(1,n);
for k = 1:length(offsetFreq)
  M = expm([0 offPhase(k) 0; -offPhase(k) 0 -flipAngle; 0 flipAngle 0])*M0;
  Mz90(k) = M(3);
end

x = offPhase/2;
sincMz = 1 - (1-cos(flipAngle))*(sin(x)./x).^2;

% plotting
%--------------------------------------------------------------------
h = plot(offsetFreq,Mz180,'b',offsetFreq,Mz90,'k',offsetFreq,sincMz,'r:');
set(h(1:2),'LineWidth',2);
xlabel('frequency offset [MHz]');
ylabel('z magnetization after pulse');
title(sprintf('rectangular pulse with length %g ns',tp*1e3));
legend('pi pulse','pi/2 pulse','sinc');
