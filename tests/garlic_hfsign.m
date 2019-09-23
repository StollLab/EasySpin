function [err,data] = test(opt,olddata)

% Test whether simulated spectrum is invariant with respect to sign
% of hyperfine coupling.

Sys.Nucs = '1H';

N = 1e5;
A = 100;

Opt.Method = 'exact';

% Frequency sweep simulation
Exp1.Field = 350;
Exp1.mwRange = [9.7 9.9];
Exp1.nPoints = N;
Sys.lwpp = 2;
Sys.A = -A; [nu,spcnu1] = garlic(Sys,Exp1,Opt);
Sys.A = +A; [nu,spcnu2] = garlic(Sys,Exp1,Opt);

% Field sweep simulation
Exp2.Range = [346 354];
Exp2.mwFreq = 9.81;
Exp2.nPoints = N;
Sys.lwpp = 0.1;
Sys.A = -A; [B,spcB1] = garlic(Sys,Exp2,Opt);
Sys.A = +A; [B,spcB2] = garlic(Sys,Exp2,Opt);

err = ~areequal(spcnu1,spcnu2,1e-3,'rel') || ~areequal(spcB1,spcB2,1e-3,'rel');

if opt.Display
  subplot(2,1,1)
  plot(B,spcnu1,B,spcnu2);
  legend('negative A','positive A');
  legend boxoff
  title('frequency sweep');
  subplot(2,1,2)
  plot(B,spcB1,B,spcB2);
  legend('negative A','positive A');
  legend boxoff
  title('field sweep');
end

data = [];
