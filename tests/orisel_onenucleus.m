function [err,data] = test(opt,olddata)

%====================================================
% Test 3: S=1/2, regression test
%====================================================
Sys = struct('S',1/2,'g',[2 2.1 2.2],'HStrain',[1 1 1]*100);
Sys = nucspinadd(Sys,'63Cu',[200 200 400]);
Exp = struct('mwFreq',9,'Field',310,'ExciteWidth',500);
Opt = struct('nKnots',91,'Display',0);
w = orisel(Sys,Exp,Opt);

data.w = w;

if ~isempty(olddata)
  err = ~areequal(olddata.w/max(olddata.w),w/max(w),1e-5,'abs');
else
  err = [];
end

if opt.Display
  w
  (olddata.w-w)./w
end