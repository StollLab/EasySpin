function ok = test()

% resfields_eig should give the same thermal intensities as resfields
% when nuclei are present (nuclear population normalization)

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = 200;  % MHz

Exp.mwFreq = 95;  % GHz
Exp.Range = [3300 3500];  % mT
Exp.Temperature = 1;  % K

Opt.Threshold = 0;

[B1,I1] = resfields(Sys,Exp,Opt);
[B2,I2] = resfields_eig(Sys,Exp,Opt);

[B1,idx] = sort(B1);
I1 = I1(idx);

ok = numel(B1)==numel(B2) && ...
  areequal(B1(:),B2(:),1e-6,'rel') && areequal(I1(:),I2(:),1e-3*max(I1),'abs');
