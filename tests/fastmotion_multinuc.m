function [err,data] = test(opt,olddata)

%======================================================
% Linewidths of hypothetical system with two 14N
%======================================================

System.g = [2.0088 2.0064 2.0027];
System.Nucs = '14N,14N';
System.A = mt2mhz([7 6 32; 5 20 4]/10); % MHz
tcorr = 1e-10; % s
field = 350; % mT

[lw,mI] = fastmotion(System,field,tcorr);

data.lw = lw;
data.mI = mI;

if ~isempty(olddata)
  ok = areequal(lw,olddata.lw,1e-6,'rel') & areequal(mI,olddata.mI);
  err = ~ok;
else
  err = [];
end

if opt.Display
  lw,mI
  olddata.lw
  olddata.mI
end
