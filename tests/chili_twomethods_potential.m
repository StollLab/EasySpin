function [err,data] = test(opt,olddata)

Sys.g = [2 2.05 2.1];
Sys.tcorr = 1e-9;
%Sys.lambda = 1;
Sys.Potential = [2 0 0 1];

Exp.Field = 350;
Exp.mwRange = [9.4 10.6];

Opt.LLKM = [6 0 2 2];
Opt.LiouvMethod = 'Freed';
[x,y] = chili(Sys,Exp,Opt);
Opt.LiouvMethod = 'general';
[x,y2] = chili(Sys,Exp,Opt);

if opt.Display
  subplot(4,1,1:3);
  scale = max(y);
  plot(x,y/scale,x,y2/scale);
  legend('Freed method','general method');
  legend boxoff
  subplot(4,1,4);
  plot(x,(y-y2)/scale);
end

err = ~areequal(y,y2,1e-5*max(y));
data = [];
