function ok = test()

% Test whether weights in global fitting work correctly.

x = linspace(0,1000);

p0 = 500;
data{1} = gaussian(x,p0-50,100);
data{2} = gaussian(x,p0+50,100);

p0 = 500;
p_lo = 100;
p_hi = 900;

FitOpt.AutoScale = 'lsq';
FitOpt.Method = 'simplex';
FitOpt.Verbosity = 0;
for i = 1:2

  FitOpt.weights = zeros(1,2);
  FitOpt.weights(i) = 1;
  result = esfit_legacy(data,@model,p0,p_lo,p_hi,FitOpt);

  ok(i) = areequal(result.fit{i}(:),data{i}(:),1e-3,'abs');

end

end

function data = model(p)

x = linspace(0,1000);
data{1} = gaussian(x,p(1),100);
data{2} = gaussian(x,p(1),100);

end
