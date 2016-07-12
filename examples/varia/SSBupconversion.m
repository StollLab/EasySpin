% Single sideband upconversion with upsampling
%==========================================================================
% In this example, a cosine wave at a low frequency is upconverted by
% mixing with an LO and the results of upconversion with a double-sideband
% is compared to that with a single-sideband mixer selecting the upper or
% lower sideband.

clear, clf

% Define input signal
%-------------------------------------------------------------
tIn = 0:0.0001:0.200; % time axis, µs
freq = 20; % input signal frequency, MHz
signalIn = cos(2*pi*freq*tIn); % input cosine function

LOfreq = 200; % LO frequency for upconversion, MHz

Opt.dtOut = 0.01e-3;
[tOut,Outboth] = rfmixer(tIn,signalIn,LOfreq,'DSB',Opt); % DSB
[tOut,Outupper] = rfmixer(tIn,signalIn,LOfreq,'SSBup',Opt); % SSB, upper sideband
[tOut,Outlower] = rfmixer(tIn,signalIn,LOfreq,'SSBdown',Opt); % SSB, lower sideband

% Plot
%--------------------------------------------------------------------
plot(tIn,signalIn,tOut,Outboth-2,tOut,Outupper-4,tOut,Outlower-6)
legend('input','upconverted - DSB','upconverted - upper SSB','upconverted - lower SSB')
xlabel('t [\mus]')
set(gca,'YTick',[])
axis tight
