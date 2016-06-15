function [err,data] = test(opt,olddata)

clear Sys Exp
Sys.S = 1;
Sys.D = -3*30e3*[1 0.01];

Sys.S = 1;
Sys.D = 30e3*rand*10*[1 0.1];
Exp.CrystalOrientation = rand(1,3)*pi;
Opt.Transitions = [1 3];

mwPolarization{1} = 'linear'; Mode{1} = {rand(1,2)*pi, [],'perpendicular','parallel'};
mwPolarization{2} = 'circular'; Mode{2} = {rand*pi,[],'perpendicular','parallel'};
mwPolarization{3} = 'circular+'; Mode{3} = {rand*pi,[],'perpendicular','parallel'};
mwPolarization{4} = 'circular-'; Mode{4} = {rand*pi,[],'perpendicular','parallel'};
mwPolarization{5} = 'unpolarized'; Mode{5} = {rand*pi,[],'perpendicular','parallel'};
mwPolarization{6} = ''; Mode{6} = {rand(1,2)*pi,[],'perpendicular','parallel'};

for k = 1:numel(mwPolarization)
  for q = 1:numel(Mode{k})
    Exp.mwPolarization = mwPolarization{k};
    Exp.Mode = Mode{k}{q};
    [a,b] = resfreqs_matrix(Sys,Exp,Opt);
  end
end

err = 0;

data = [];
