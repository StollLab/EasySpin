function [err,data] = test(opt,olddata)
rand('twister',5);
Sys.S = randi_(6)/2;
Sys.g = rand(1,3)*4;
Sys.D = rand(1,2)* 10*30e3;

Exp.Temperature = [1:5:300];
Exp.Field = 0:100:500;
nFields = length(Exp.Field);

keywords = {'MvsB', 'MvsBCGS', 'MvsBSI' ...
  'Chi', 'ChiT', '1overChi', 'MuEff', ...
  'ChiSI', 'ChiTSI', '1overChiSI', 'MuEffSI', ...
  'ChiCGS', 'ChiTCGS', '1overChiCGS', 'MuEffCGS'};

mu1 = curry(Sys,Exp);
[mu,chi] = curry(Sys,Exp);

fullerr = @(a,b)abs(sum(a(:)-b(:)));
err = fullerr(mu1,mu(:));

for k = 1:15
  Opt.Output = keywords{k};
  out{k} = curry(Sys,Exp,Opt);
  
  switch k
    case {1,2}, test = mu; %MvsB,MvsBCGS
    case 3, test = mu *avogadro* bmagn; %MvsBSi
    case {4,8},test = chi; %Chi,ChiSI
    case {5,9}, test = chi.*repmat(Exp.Temperature(:).',nFields,1); %ChiT,ChiTSi
    case {6,10}, test = 1./chi; %1overChi,1overChiSI
    case {7,11}, test = sqrt(chi.*repmat(Exp.Temperature(:).',nFields,1)*8); %MuEff, MuEffSI
    case 12, test = chi/(4*pi*1e-6); %ChiCGS
    case 13, test = chi.*repmat(Exp.Temperature(:).',nFields,1)/(4*pi*1e-6); %ChiTCGS
    case 14, test = 1./chi*(4*pi*1e-6); %1overChiCGS
    case 15, test = sqrt(chi.*repmat(Exp.Temperature(:).',nFields,1)*8/(4*pi*1e-6)); %MuEffCGS
  end
  err = err + fullerr(out{k},test);
end


for l = 1:10
  r = randi_(10);
  Opt.Output = '';
  for x = 1:r
    t(x) = randi_(15);
    Opt.Output = [Opt.Output,' ', keywords{t(x)}];
  end
  [test2{1:r}] = curry(Sys,Exp,Opt);
  for x = 1:r
    err = err + fullerr(out{t(x)},test2{x});
  end
end

err = err > 1e-14;

data =[];


