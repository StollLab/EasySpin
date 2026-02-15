function ok = test(opt)

% Test whether tilting a tensor and the integration grid leaves the
% simulated spectrum invariant

Sys.g = [1.9 2 2.1];
Sys.gFrame = [0 0 0];

Exp.mwFreq = 9.5;

spc0 = pepper(Sys,Exp);

Sys.gFrame = deg2rad([10 20 30]);

Opt.GridSymmetry = 'D2h';
Opt.GridFrame = Sys.gFrame;

[B,spc1] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(B,spc0,B,spc1);
  legend('untiled gFrame, automatic grid','tilted gFrame, manually tilted grid');
end

ok = areequal(spc0,spc1,1e-5,'rel');

end
