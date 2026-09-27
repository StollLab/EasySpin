function ok = test()

% Without temperature, polarizations in frequency sweeps are proportional
% to the transition frequency, as in the high-temperature limit

Sys.S = 1;
Sys.D = 3000;  % MHz

Exp.Field = 350;  % mT
Exp.SampleFrame = [0 0 0];

[nu,Int0] = resfreqs_matrix(Sys,Exp);
Exp.Temperature = 1e4;  % K
[~,Int1] = resfreqs_matrix(Sys,Exp);

ok(1) = numel(nu)==2 && abs(diff(nu))>1e3;
ok(2) = areequal(Int0(1)/Int0(2),Int1(1)/Int1(2),1e-3,'rel');
