function ok = test()

% Test whether FitOpt.AutoScale and FitOpt.BaseLine work correctly.

p = [400 600 300];
[x,data0] = model(p);


for k = 1:numel(data0)
  noiselevel = max(data0{k})*0.1;
  offset = rand*max(data0{k});
  noise{k} = randn(size(data0{k}))*noiselevel;
  data0{k} = data0{k} + noise{k} + offset;
end

p0 = [500 100 100];
plo = [0 10 10];
pup = [1e3 500 500];

autoscaleopts = {'none' 'lsq'};
baselineopts1 = {0 1 2};
baselineopts2 = {0 1 2};

for k = 1:numel(baselineopts1)

  for j = 1:numel(baselineopts2)

    FitOpt.BaseLine = [baselineopts1{k} baselineopts2{j}];

    for ii = 1:numel(data0)
      switch FitOpt.BaseLine(ii)
        case 0
          bl = 0.1;
        case 1
          bl = 0.1 + rand*1e-4*x;
        case 2
          bl = 0.1 + rand*1e-5*(x-300) + rand*1e-7*x.^2;
      end
      data{ii} = data0{ii}+bl;
    end

    for i = 1:numel(autoscaleopts)

      FitOpt.Verbosity = 0;
      FitOpt.AutoScale = autoscaleopts{i};
      result = esfit_legacy(data,@model,p0,plo,pup,FitOpt);

      noisevec = [noise{1} noise{2}];
      ok(i,j,k) = areequal(std(noisevec),result.rmsd,0.2,'rel');

    end
  end
end

end

function [x,data] = model(p)

x = linspace(0,1e3,1e3);

data{1} = gaussian(x,p(1),p(2));
data{2} = gaussian(x,p(1),p(3));

end
