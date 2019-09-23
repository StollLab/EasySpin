function [err,data] = test(opt,olddata)

%=======================================================
% Test 2: S=1/2 + 1H + 14N
%=======================================================

Sys = struct('S',1/2,'g',[2 2 2.3],'Nucs','14N,1H','A',[1 1 1;1 1 1]*20);
Exp = struct('mwFreq',9.5,'Range',[200 400]);
Exp.CrystalOrientation = [10 10 0]*pi/180;
Opt = struct('Threshold',0,'Verbosity',0,'Perturb',0);
Opt.RejectionRatio = 1e-3;

% should give 36 transitions
[B,A] = eigfields(Sys,Exp,Opt);

data.B = B;
data.A = A;

if (opt.Display)
  plot(B,A,'o',olddata.B,olddata.A,'.');
  legend('new','old');
end

if ~isempty(olddata)
  thr = 1e-4;
  err = ~areequal(olddata.B,B,thr,'rel') | ~areequal(olddata.A,A,thr,'rel');
else
  err = [];
end

