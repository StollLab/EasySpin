function ok = test()

% Add multiple nuclei

Sys.S = 1/2;
Sys.g = [2.0054 2.0042 2.0022];
Sys = nucspinadd(Sys,'1H',[1 2 3]);
Sys = nucspinadd(Sys,'1H',[4 5 6]);
Sys = nucspinadd(Sys,'17O',[7 8 9],[],[-1 -1 2],[]);

ok = numel(Sys.Q)==9;
