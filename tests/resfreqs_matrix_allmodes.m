function ok = test()

Sys.S = 1;
Sys.D = -3*30e3*[1 0.01];

Sys.S = 1;
Sys.D = 30e3*rand*10*[1 0.1];
Exp.CrystalOrientation = rand(1,3)*pi;
Opt.Transitions = [1 3];

Mode{1} = 'perpendicular';
Mode{2} = 'parallel';
Mode{3} = {pi/6, pi/4};
Mode{4} = {pi/6, 'circular+'};
Mode{5} = {pi/6, 'circular-'};
Mode{6} = {pi/6, 'unpolarized'};

for k = 1:numel(Mode)
  [~,~] = resfreqs_matrix(Sys,Exp,Opt);
end

ok = true;
