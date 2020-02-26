function ok = test(opt)

% Pseudo-Voigtian

x = linspace(-200,200,1e4); x0 = 34;
wG = 20; wL = 15;
dif = 0;
alphaGauss = 0.4;

yG =   gaussian(x,x0,wG,dif);
yL = lorentzian(x,x0,wL,dif);

yA = lshape(x,x0,[wG wL],dif,alphaGauss);
yB = yG*alphaGauss + yL*(1-alphaGauss);

if (opt.Display)
  plot(x,yA,'b',x,yB,'r');
  legend('lshape','manual sum');
end

ok = areequal(yA,yB,1e-10,'rel');
