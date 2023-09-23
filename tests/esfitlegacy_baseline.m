function ok = test()

% Test whether FitOpt.AutoScale and FitOpt.BaseLine work correctly.

x = linspace(300,400,100);
pos = 350;
lw = 7;
model = @(p) 5*exp(-(x-p(1)).^2/(2*p(2)^2));

spcsim = model([pos lw]);
noiselevel = max(spcsim)*0.1;
noisevec = randn(size(spcsim))*noiselevel;
spcexp0 = spcsim + noisevec;

p0 = [340 5];
pv = [50 20];

autoscaleopts = {'none' 'lsq' 'maxabs'};
baselineopts = {[] 0 1 2};

for j = 1:numel(baselineopts)
  
  FitOpt.BaseLine = baselineopts{j};
  
  if isempty(FitOpt.BaseLine)
    bl = 0;
  else
    switch FitOpt.BaseLine
      case 0
        bl = 0.1;
      case 1
        bl = 0.1 + 0.01*(x-300);
      case 2
        bl = 0.1 + 0.01*(x-300) + 0.0002*(x-340).^2;
    end
  end
  spcexp = spcexp0+bl;
  
  for i = 1:numel(autoscaleopts)
    
    FitOpt.Verbosity = 0;
    FitOpt.AutoScale = autoscaleopts{i};
    result = esfit_legacy(spcexp,model,p0,pv,FitOpt);
    
    ok(i,j) = areequal(std(noisevec),result.rmsd,0.2,'rel');
    
  end
end

end

