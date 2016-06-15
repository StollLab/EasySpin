function [err,data] = resfields_sparsefull(opt,olddata)

% sparse vs full matrics

Sys.S = [3/2 1];
Sys.D = 1000*[1 1.01];
Sys.ee = 202;

Exp.Field = 10000;

Opt.SparseMode = true;
[B1,I1] = resfreqs_matrix(Sys,Exp,Opt);
data1 = [B1, I1];

Opt.SparseMode = false;
[B2,I2] = resfreqs_matrix(Sys,Exp,Opt);
data2 = [B2,I2];

err = any(abs(data1(:)-data2(:))>1e-2);
data = [];
