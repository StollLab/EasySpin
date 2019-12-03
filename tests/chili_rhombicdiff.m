function [err,data] = test(opt,olddata)

%=======================================================
% Rhombic diffusion tensor
%=======================================================

Sys = struct('g',[2.1 2.0 1.9],'lw',1);
Exp = struct('mwFreq',9.5,'CenterSweep',[340 100]);

q = 0.3e8;
Sys.Diff = [10  1  1]*q; [x,y(1,:)] = chili(Sys,Exp);
Sys.Diff = [1 10 1]*q; [x,y(2,:)] = chili(Sys,Exp);
Sys.Diff = [1 1  10]*q; [x,y(3,:)] = chili(Sys,Exp);

for k=1:3
  y(k,:) = y(k,:)/max(y(k,:));
end

data.x = x;
data.y = y;

if ~isempty(olddata)
  if opt.Display
    plot(data.x,data.y,'r',data.x,olddata.y,'g');
  end
  ok = areequal(y,olddata.y,1e-3,'abs');
  err = ~ok;
else
  err = [];
end
