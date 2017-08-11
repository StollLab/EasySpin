function [err,data] = test(opt,olddata)

rand('twister',5);

Sys.S = randi_(6)/2;
Sys.g = rand(1,3)*4;
Sys.D = rand(1,2)* 10*30e3;

T = 1:10:300;
Exp.Temperature = T;
Exp.Field = 0:100:500;
nFields = length(Exp.Field);

keywords = {'mu', 'mumol', 'muBM', 'chimol', 'chimolT', '1/chimol', 'mueff'};

% Calculate basic quantities in SI units
Opt.Output = 'mu chimol';
Opt.Units = 'SI';
[mu_SI, chimol_SI] = curry(Sys,Exp,Opt);

% Make some basic conversions
% (see Hatscher et al, Pure Appl. Chem. 77(2), 497-511 (2005), Table 1
mu_CGS = mu_SI/1e-3;
chimol_CGS = chimol_SI/(4*pi*1e-6);

chimolT_SI = chimol_SI.*repmat(T(:).',nFields,1);
chimolT_CGS = chimol_CGS.*repmat(T(:).',nFields,1);

muB_SI = bmagn;
muB_CGS = muB_SI/1e-3;
muBM_SI = mu_SI/muB_SI;
muBM_CGS = mu_CGS/muB_CGS;

fullerr = @(a,b)abs(max(a(:)-b(:)));
err = 0;

for k = 1:numel(keywords)
  
  Opt.Output = keywords{k};
  
  Opt.Units = 'SI';
  out = curry(Sys,Exp,Opt);
  switch keywords{k}
    case 'mu', test = mu_SI;
    case 'muBM', test = muBM_SI;
    case 'mumol', test = mu_SI*avogadro;
    case 'chimol', test = chimol_SI;
    case 'chimolT', test = chimolT_SI;
    case '1/chimol', test = 1./chimol_SI;
    case 'mueff', c = 3*boltzm/avogadro/bmagn^2/mu0; test = sqrt(chimolT_SI*c);
  end
  err = err + fullerr(out,test);
  
  Opt.Units = 'CGS';
  out = curry(Sys,Exp,Opt);
  switch keywords{k}
    case 'mu', test = mu_CGS;
    case 'muBM', test = muBM_CGS;
    case 'mumol', test = mu_CGS*avogadro;
    case 'chimol', test = chimol_CGS;
    case 'chimolT', test = chimolT_CGS;
    case '1/chimol', test = 1./chimol_CGS;
    case 'mueff', c = 1/0.1250494086; test = sqrt(chimolT_CGS*c);
  end
  err = err + fullerr(out,test);
  
end

err = err > 1e-12;

data = [];
