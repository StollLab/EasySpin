function ok = test(opt)

% Comparison of axial and orthorhombic spectra
% with interpolation on and off (same total number
% of knots).

Sys = struct('S',.5,'g',[1.9,2.01,2.3],'lw',1);
Exp = struct('Range',[285 365],'mwFreq',9.5);
Opt = struct('Verbosity',opt.Verbosity);

n1 = 19; k1 = 5;
n2 = (n1-1)*k1 + 1; k2 = 1;
Opt.GridSymmetry ='D2h';
Opt.GridSize = [n1 k1]; [x,y1] = pepper(Sys,Exp,Opt);
Opt.GridSize = [n2 k2]; [x,y2] = pepper(Sys,Exp,Opt);
Opt.GridSymmetry = 'Dinfh';
Opt.GridSize = [n1 k1]; [x,y3] = pepper(Sys,Exp,Opt);
Opt.GridSize = [n2 k2]; [x,y4] = pepper(Sys,Exp,Opt);

ok = true;

if (opt.Display)
  title('Test 7: Interpolation on/off');
  subplot(2,1,1); plot(x,y1/max(y1),x,y2/max(y2));
  legend('on','off','Location','BestOutside');
  subplot(2,1,2); plot(x,y3/max(y3),x,y4/max(y4));
  legend('on','off','Location','BestOutSide');
end
