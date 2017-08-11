function [err,data] = test(opt,olddata)

rand('twister',5);

Sys.S = randi_(6)/2;
Sys.g = rand(1,3)*4;
Sys.D = rand(1,2)* 10*30e3;


keywords = {'mu', 'chimol', 'chimolT', '1/chimol', 'mueff'};

err = 0;

for k = 1:numel(keywords)
  Opt.Output = ['onecolumn ' keywords{k}];
  Exp.Temperature = [1:10];
  Exp.Field = 100 * ones(5,1);
  nTemps = length(Exp.Temperature);
  nFields = length(Exp.Field);
  t2 = curry(Sys,Exp,Opt);
  
  
  for m = 1:(nTemps-1)
    for n = 2:nFields
      err = err +abs(t2((m-1)*nFields+1)-t2((m-1)*nFields+n));
    end
  end
  
  Exp.Temperature = 100 * ones(5,1);
  Exp.Field = 100 * 1:10;
  nTemps = length(Exp.Temperature);
  nFields = length(Exp.Field);
  t2 = curry(Sys,Exp,Opt);
  
  for n = 1:nFields
    for m = 2:(nTemps-1)
      err = err +abs(t2((n-1)*nTemps+1)-t2((n-1)*nTemps+m));
    end
  end
end

data = [];


