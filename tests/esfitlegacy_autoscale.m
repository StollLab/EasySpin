function ok = test()

% Test whether FitOpt.AutoScale works correctly.

t = linspace(0,20,501);
modelfcn = @(p) exp(-p(1)*t)+p(2);
pmodel = [0.3 0.2];
scale = 5.34564564;
yexp = scale*modelfcn(pmodel);

p0 = [0.1 0.4];
lb = [0 0];
ub = [1 1];

FitOpt = struct;
FitOpt.Method = 'simplex fcn';
FitOpt.Verbosity = 0;
FitOpt.AutoScale = 'lsq';

result = esfit_legacy(yexp,modelfcn,p0,lb,ub,FitOpt);

ok(1) = areequal(result.pfit(:),pmodel(:),1e-3,'abs');
ok(2) = areequal(result.scale,scale,1e-3,'abs');

end
