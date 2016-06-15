function [err,data] = test(opt,olddata)

% Test transition selection for resfreqs_matrix

Sys.S = 1;
Sys.D = 3*30e3*[1 0.2];
Sys.lwpp = 2000;

Exp.Range = [0 3000];
nu = resfreqs_matrix(Sys,Exp);
err(1) = numel(nu)~=3;

Exp.Range = [30 80];
nu = resfreqs_matrix(Sys,Exp);
err(2) = numel(nu)~=2;

Exp.Range = [100 200];
nu = resfreqs_matrix(Sys,Exp);
err(3) = numel(nu)~=1;

Exp.Range = [200 400];
nu = resfreqs_matrix(Sys,Exp);
err(4) = numel(nu)~=0;

Exp.Range = [10 20];
nu = resfreqs_matrix(Sys,Exp);
err(5) = numel(nu)~=0;

err = any(err);
data = [];
