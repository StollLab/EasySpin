function ok = test()

% Test whether chili returns the correct number of subspectra
% when Opt.Output='components'.

clear, clc

Sys1.g = [2 2.1 2.2];
Sys1.tcorr = 1e-9;  % s

Sys2.g = [2 2.1 2.2];
Sys2.tcorr = 0.3e-9;  % s

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [280 360];  % mT

Opt.Output = 'components';

[B,spc] = chili(Sys1,Exp,Opt);
ok(1) = size(spc,1)==1;

[B,spc] = chili(Sys2,Exp,Opt);
ok(2) = size(spc,1)==1;

[B,spc] = chili({Sys1,Sys2},Exp,Opt);
ok(3) = size(spc,1)==2;

end
