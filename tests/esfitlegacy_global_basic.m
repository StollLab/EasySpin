function ok = test(opt)

% Check that global fitting with multiple data input

p = [400 100 10];
[x,modeldata] = model(p);

p0 = [500 50 50];
plo = [100 1 1];
pup = [900 200 200];

FitOpt.Verbosity = 0;
result = esfit_legacy(modeldata,@model,p0,plo,pup,FitOpt);

ok = result.rmsd/max([result.fit{1}; result.fit{2}])<3e-2;

if opt.Display
  clf
  for i = 1:numel(modeldata)
    nexttile
    plot(x,modeldata{i},x,result.fit{i});
    legend('exp','fit');
  end
end

end

function [x,data] = model(p)

x = linspace(0,1e3,1e3);

data{1} = gaussian(x,p(1),p(2));
data{2} = gaussian(x,p(1),p(3));

end