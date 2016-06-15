function [err,data] = test(opt,olddata)

%======================================================
% Linewdiths of hypothetical system with two 14N
%======================================================

System = struct('g',[2.0088 2.0064 2.0027],'Nucs','14N,14N');
System.A = mt2mhz([7 6 32; 5 20 4]/10);
tcorr = 1e-10;
field = 350;

[lw,mI] = fastmotion(System,field,tcorr);

data.lw = lw;
data.mI = mI;

if ~isempty(olddata)
  ok = areequal(lw,olddata.lw,1e-6) & areequal(mI,olddata.mI);
  err = ~ok;
else
  err = [];
end

if opt.Display
  lw,mI
  olddata.lw
  olddata.mI
end
