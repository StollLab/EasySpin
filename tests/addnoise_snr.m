function ok = test()

% Make sure the SNR is as defined

signal_amplitude = 123;
SNR = 20;

y = zeros(1,1e4);
y(1) = signal_amplitude;

yf = addnoise(y,SNR,'f'); yf(1) = [];
yu = addnoise(y,SNR,'u'); yu(1) = [];
yn = addnoise(y,SNR,'n'); yn(1) = [];

noise_stddev = [std(yf) std(yu) std(yn)];
SNReff = signal_amplitude./noise_stddev;

ok = areequal(SNReff,SNR*[1 1 1],0.05,'rel');
