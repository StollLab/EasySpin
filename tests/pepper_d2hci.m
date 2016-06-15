function [err,data] = test(opt,olddata)

% D2h vs Ci computation

Sys = struct('S',.5,'g',[1.9,2,2.3],'Nucs','1H','A',[20 80 ...
      300],'AFrame',[-pi/4 -pi/2 0]);
Pg = symm(Sys);
%fprintf('Correct symmetry is %s\n',Pg);
Exp = struct('Range',[285 365],'mwFreq',9.5,'nPoints',2e3,'Harmonic',0);
Opt = struct('nKnots',[20 4],'Method','perturb');
Symmetry = {'Ci','C2h','D2h'};

for k = 1:length(Symmetry)
  Opt.Symmetry = Symmetry{k};
  [x,y(k,:)] = pepper(Sys,Exp,Opt);
end

err = 5e-3*max(y(1,:));
err = ~areequal(y(1,:),y(2,:),err);

data = [];

if (opt.Display)
  title('D2h, C2h, Ci computation (D2h is incorrect)');
  plot(x,y);
  legend(Symmetry{:});
end
