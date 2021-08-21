% ESEEM dead time correction with ctafft()
%==========================================================================
clear, clf

% Construct a time-domain signal containing some frequency components
dt = 0.01; tmax = 3; t = 0:dt:tmax;
Freq = [22 35 40 56 76];
Ampl = [0.1 1 1 1 1];
decay = [4 1 0.7 1 0.5]*0.05;
td = zeros(1,length(t));
td0 = 0;
for i = 1:length(Freq)
  td0 = td0 + Ampl(i)*exp(2i*pi*Freq(i)*t).*exp(-t/decay(i));
end

% Various Fourier transforms of the signal
fd0 = fft(td0,1024);
td_d = td0(3:end);
fd_d = fft(td_d,1024);
fd_cta = ctafft(td_d,40,1024);

%Graphical rendering
subplot(3,1,1); plot(abs(fd0),'k'); axis tight;
title('magnitude, no dead time');
subplot(3,1,2); plot(abs(fd_d),'r'); axis tight;
title('magnitude, with dead time');
subplot(3,1,3); plot(fd_cta); axis tight;
title('cross term averaged, with dead time');
