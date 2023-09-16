function ok = test(opt)

x = linspace(0,100);
data = gaussian(x,30,10);

baseline0 = 0.1+0.002*x - (0.002*x).^2;
data = data + baseline0;

region = x<10 | x>50;
[~,b] = basecorr(data,2,2,region);

if opt.Display
  plot(x,data,x,baseline0,x,b,'--');
  xregion(0,10);
  xregion(50,100);
  legend('data','synthetic baseline','fitted baseline');
end

ok = areequal(b,baseline0,1e-7,'abs');
