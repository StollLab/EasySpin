function [err,data] = test(opt,olddata)

% adding nuclei with full hyperfine matrices
%====================================================
Afull = rand(3);

% no nucleus present
Sys1.g = 2;
Sys1 = nucspinadd(Sys1,'14N',Afull);
err(1) = (numel(Sys1.A)~=9);

% one nucleus present, isotropic hfi
Sys2.Nucs = '1H'; Sys2.A = 1;
Sys2 = nucspinadd(Sys2,'14N',Afull);
err(2) = (numel(Sys2.A)~=18);

% one nucleus present, axial hfi
Sys3.Nucs = '1H'; Sys3.A = [1 2];
Sys3 = nucspinadd(Sys3,'14N',Afull);
err(3) = (numel(Sys3.A)~=18);

% one nucleus present, rhombic hfi
Sys4.Nucs = '1H'; Sys4.A = [1 2 3];
Sys4 = nucspinadd(Sys4,'14N',Afull);
err(4) = (numel(Sys4.A)~=18);

% one nucleus present, full hfi
Sys5.Nucs = '1H'; Sys5.A = [1 2 3; 4 5 6; 7 8 9];
Sys5 = nucspinadd(Sys5,'14N',Afull);
err(5) = (numel(Sys5.A)~=18);

if any(err)
  err
end

err = any(err);
data = [];
