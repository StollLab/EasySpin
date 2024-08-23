function ok = test(opt)

% Compare Voigtian dispersion calculated two ways

x = linspace(300,400,1001);
x0 = 340;

% Case 1: Gaussian FWHM < Lorentzian FWHM
fwhmGL = [1 2];
ydisp1a = voigtian(x,x0,fwhmGL,0,pi/2);
[~, ydisp2a] = voigtian(x,x0,fwhmGL,0,0);

% Case 2: Gaussian FWHM > Lorentzian FWHM
fwhmGL = [1 0.7];
ydisp1b = voigtian(x,x0,fwhmGL,0,pi/2);
[~, ydisp2b] = voigtian(x,x0,fwhmGL,0,0);

% Plotting
if opt.Display
  subplot(2,1,1)
  plot(x,ydisp1a,x,ydisp2a,'.');
  legend('direct','two output args');
  title('fwhmG<fwhmL')
  subplot(2,1,2)
  plot(x,ydisp1b,x,ydisp2b,'.');
  legend('direct','two output args');
  title('fwhmG>fwhmL')
end

threshold = 1e-8;
ok(1) = areequal(ydisp1a,ydisp2a,threshold,'rel');
ok(2) = areequal(ydisp1b,ydisp2b,threshold,'rel');
