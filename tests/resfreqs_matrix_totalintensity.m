function [err,data] = test(opt,olddata)

% Total transition intensity should be independent of orientation of k
% vector relative to B0, no matter what the crystal orientation and the
% polarization type and angle are.

clear Sys Exp
Sys.S = 1;
Sys.D = 30e3*rand*10*[1 1/3];
Exp.CrystalOrientation = rand(1,3)*pi;
Opt.Transitions = [ 1 2; 1 3; 2 3];

for n=10:-1:1
  Exp.mwPolarization = 'linear';
  Exp.Mode = [pi/2 rand*pi];
  [p(n,:),i(n,:)] = resfreqs_matrix(Sys,Exp,Opt);
end
err(1) = any(abs(sum(i')-sum(i(1,:)))>1e-10);

for n=10:-1:1
  Exp.mwPolarization = 'linear';
  Exp.Mode = [rand*pi pi/2];
  [p2(n,:),i2(n,:)] = resfreqs_matrix(Sys,Exp,Opt);
end
err(2)= any(abs(sum(i2')-sum(i2(1,:)))>1e-10);

for n=10:-1:1
  Exp.mwPolarization = 'unpolarized';
  Exp.Mode =  rand*pi;
  [p3(n,:),i3(n,:)] = resfreqs_matrix(Sys,Exp,Opt);
end
err(3) = any(abs(sum(i3')-sum(i3(1,:)))>1e-10);

for n=10:-1:1
  Exp.mwPolarization = 'circular-';
  Exp.Mode =  rand*pi;
  [p4(n,:),i4(n,:)] = resfreqs_matrix(Sys,Exp,Opt);
end
err(4) = any(abs(sum(i4')-sum(i4(1,:)))>1e-10);

for n=10:-1:1
  Exp.mwPolarization = 'circular+';
  Exp.Mode =  rand*pi;
  [p5(n,:),i5(n,:)] = resfreqs_matrix(Sys,Exp,Opt);
end
err(5) = any(abs(sum(i5')-sum(i5(1,:)))>1e-10);

err = any(err);

data = [];
