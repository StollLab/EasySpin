function ok = test()

% Test whether garlic returns the correct number of subspectra
% when Opt.Output='components'.

Sys1.g = 2;
Sys1.Nucs = '1H';
Sys1.A = 10;
Sys1.n = 3;
Sys1.lwpp = 0.03;

Sys2.g = 2;
Sys2.Nucs = 'H';
Sys2.A = 10;
Sys2.n = 3;
Sys2.lwpp = 0.03;

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [336 343];  % mT

Opt.Output = 'components';

[B,spc] = garlic(Sys1,Exp,Opt);
ok(1) = size(spc,1)==1;

[B,spc] = garlic(Sys2,Exp,Opt);
ok(2) = size(spc,1)==2;

[B,spc] = garlic({Sys1,Sys2},Exp,Opt);
ok(2) = size(spc,1)==3;

end
