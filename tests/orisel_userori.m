function [err,data] = test(opt,olddata)

%====================================================
% Test 2: user-supplied orientations
%====================================================
Sys = struct('S',1/2,'g',[2 2 2.2]);
Sys.HStrain = [1 1 1]*50;
Exp = struct('mwFreq',9.5,'Field',320);
Exp.ExciteWidth = 100;
N = 265;
Exp.Orientations = [zeros(1,N); linspace(0,pi/2,N)];
w = orisel(Sys,Exp);
err = (numel(w)~=N);

data = [];
