function ok = test(opt)

sys.Nucs = '1H,1H';
sys.A = [3 0.5];
sys.lwEndor = 0.1;

exp.Sequence = 'MimsENDOR';
exp.tau = 0.01;
exp.Field = 350;
exp.Range = [.5 18];
exp.CrystalOrientation = [0 0 0];
exp.nPoints = 1024;

opt.Verbosity = 0;

opt.EndorMethod = 1;
opt.ProductRule = 0;
[x,y0] = saffron(sys,exp,opt);
opt.ProductRule = 1;
[x,y1] = saffron(sys,exp,opt);

n = @(y)y/max(abs(y));
y0 = n(y0);
y1 = n(y1);

if (opt.Display)
  plot(x,y1,x,y0,'r');
  xlabel('Frequency /MHz')
end

ok = sum(abs(y0-y1)) < 1e-12;
