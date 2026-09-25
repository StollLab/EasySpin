function ok = test()

% resfields_eig should give the same thermal intensities as resfields
% (S=1 with large zero-field splitting at low temperature)

Sys.S = 1;
Sys.g = 2;
Sys.D = [3000 500];  % MHz

Exp.mwFreq = 95;  % GHz
Exp.Range = [0 6000];  % mT
Exp.SampleFrame = [0 pi/5 pi/7];
Exp.Temperature = 2;  % K

Opt.Threshold = 0;

[B1,I1] = resfields(Sys,Exp,Opt);
[B2,I2] = resfields_eig(Sys,Exp,Opt);

[B1,idx] = sort(B1);
I1 = I1(idx);

ok = numel(B1)==numel(B2) && ...
  areequal(B1(:),B2(:),1e-6,'rel') && areequal(I1(:),I2(:),1e-3*max(I1),'abs');
