function ok = test()

% sparse vs full matrics

Sys.S = [3/2 1];
Sys.D = 1000*[1 1.01];
Sys.ee = 202;

Exp.mwFreq = 9.5;
Exp.Range = [20 1000];

Opt.SparseMode = true;
[B1,I1] = resfields(Sys,Exp,Opt);
data1 = [B1, I1];

Opt.SparseMode = false;
[B2,I2] = resfields(Sys,Exp,Opt);
data2 = [B2,I2];

ok = all(abs(data1(:)-data2(:))<1e-2);
