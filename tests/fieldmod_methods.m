function ok = test(opt)

% Assert that the fft and the conv method give similar results.

B0 = 340;  % mT
B = 300:0.01:400;  % mT
fwhm = 1;  % mT

modAmp = 20;  % mT
harmonic = 1;

spc = gaussian(B,B0,fwhm);

spc_fft = fieldmod(B,spc,modAmp,harmonic,'fft');
spc_conv = fieldmod(B,spc,modAmp,harmonic,'conv');

ok = areequal(max(spc_fft),max(spc_conv),0.02,'rel');

if opt.Display
  plot(B,spc_fft,B,spc_conv,'.');
  legend('fft','conv');
  legend boxoff
  xlabel('field (mT)');
end
