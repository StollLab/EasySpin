function ok = test()

% Make sure the SNR is as defined

signal_amplitude = 123*(1+0.1i);
SNR = 20;

y = zeros(1,1e4);
y(1) = signal_amplitude;

yf = addnoise(y,SNR,'f'); yf(1) = [];
yu = addnoise(y,SNR,'u'); yu(1) = [];
yn = addnoise(y,SNR,'n'); yn(1) = [];

noise_stddev_re = [std(real(yf)) std(real(yu)) std(real(yn))];
noise_stddev_im = [std(imag(yf)) std(imag(yu)) std(imag(yn))];

SNReff_re = real(signal_amplitude)./noise_stddev_re;
SNReff_im = imag(signal_amplitude)./noise_stddev_im;
SNReff = max([SNReff_re SNReff_im]);

ok = areequal(SNReff,SNR,0.05,'rel');
