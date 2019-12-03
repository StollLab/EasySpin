function [err,data] = test(opt,olddata)

% adding nuclei with full quadrupole matrices
%====================================================
Afull = rand(3);
Qfull = rand(3);

% no nucleus
%----------------------------------------------------
Sys1.g = 2;
Sys1 = nucspinadd(Sys1,'14N',Afull,[],Qfull,[]);
err(1) = numel(Sys1.Q)~=9;

% one nucleus, no  nqi
%----------------------------------------------------
% isotropic hfi
Sys2.Nucs = '1H'; Sys2.A = 1;
Sys2 = nucspinadd(Sys2,'14N',Afull,[],Qfull,[]);
err(2) = (numel(Sys2.Q)~=18);

% axial hfi
Sys3.Nucs = '1H'; Sys3.A = [1 2];
Sys3 = nucspinadd(Sys3,'14N',Afull,[],Qfull,[]);
err(3) = (numel(Sys3.Q)~=18);

% rhombic hfi
Sys4.Nucs = '1H'; Sys4.A = [1 2 3];
Sys4 = nucspinadd(Sys4,'14N',Afull,[],Qfull,[]);
err(4) = (numel(Sys4.Q)~=18);

% full hfi
Sys5.Nucs = '1H'; Sys5.A = [1 2 3; 4 5 6; 7 8 9];
Sys5 = nucspinadd(Sys5,'14N',Afull,[],Qfull,[]);
err(5) = (numel(Sys5.Q)~=18);

% one nucleus, with nqi principal values
%----------------------------------------------------
clear Sys*
% isotropic hfi
Sys2.Nucs = '2H'; Sys2.A = 1; Sys2.Q = [1 2 3];
Sys2 = nucspinadd(Sys2,'14N',Afull,[],Qfull,[]);
err(2) = (numel(Sys2.Q)~=18);

% axial hfi
Sys3.Nucs = '1H'; Sys3.A = [1 2];  Sys3.Q = [1 2 3];
Sys3 = nucspinadd(Sys3,'14N',Afull,[],Qfull,[]);
err(3) = (numel(Sys3.Q)~=18);

% rhombic hfi
Sys4.Nucs = '1H'; Sys4.A = [1 2 3];  Sys4.Q = [1 2 3];
Sys4 = nucspinadd(Sys4,'14N',Afull,[],Qfull,[]);
err(4) = (numel(Sys4.Q)~=18);

% full hfi
Sys5.Nucs = '1H'; Sys5.A = [1 2 3; 4 5 6; 7 8 9];  Sys5.Q = [1 2 3];
Sys5 = nucspinadd(Sys5,'14N',Afull,[],Qfull,[]);
err(5) = (numel(Sys5.Q)~=18);

if opt.Display
  err
end

err = any(err);
data = [];
