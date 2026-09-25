function ok = test()

% Separate transition spectra over a field range across zero: resonances of
% the same transition at positive and negative fields are combined

Sys.S = 1;
Sys.D = 300;
Sys.lw = 5;
Exp.mwFreq = 9.5;
Exp.Range = [-500 500];
Opt.separate = 'transitions';

[~,ysep,info] = pepper(Sys,Exp,Opt);
[~,ysum] = pepper(Sys,Exp);

nTransitions = size(unique(info.Transitions,'rows'),1);

ok = size(ysep,1)==nTransitions && areequal(sum(ysep,1),ysum,1e-8,'rel');
