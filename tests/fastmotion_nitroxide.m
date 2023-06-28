function [ok,data] = test(opt,olddata)

%======================================================
% Linewdiths of a nitroxide radical
%======================================================

System = struct('g',[2.0088 2.0064 2.0027],'Nucs','14N');
System.A = unitconvert([7.59 5.95 31.76]/10,'mT->MHz');
tcorr = 1e-10;
field = 350;

[lw,mI] = fastmotion(System,field,tcorr);

data.lw = lw;
data.mI = mI;

if ~isempty(olddata)
  if opt.Display
    lw
    olddata.lw
  end
  ok = areequal(lw,olddata.lw,1e-4,'rel') & areequal(mI,olddata.mI);
else
  ok = [];
end
