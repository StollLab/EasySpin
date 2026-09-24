function ok = test()

% Read ORCA 6 output from a parameter scan (%paras) of the hydroxyl
% radical: single-point calculations at O-H distances from 1.0 to 1.3
% Angstrom in 4 steps, with EPR properties calculated at each step.
%   ok(1)  main output file (.out)
%   ok(2)  property file (.property.txt)
%   ok(3)  main output and property file agree for each step

name = 'orca/v6.1.0/hydroxyl_parascan';
rOH = [1.0 1.1 1.2 1.3];  % Angstrom

SysM = tryload([name '.out']);
SysP = tryload([name '.property.txt']);

ok(1) = checkscan(SysM,rOH);
ok(2) = checkscan(SysP,rOH);
ok(3) = ok(1) && ok(2) && comparescans(SysM,SysP);

end

%-------------------------------------------------------------------------------
% Load file, return [] if reading fails
function Sys = tryload(fileName)
try
  Sys = orca2easyspin(fileName);
catch
  Sys = [];
end
end

%-------------------------------------------------------------------------------
% Check spin system array against the scan specification
function ok = checkscan(Sys,rOH)
nSteps = numel(rOH);
ok = numel(Sys)==nSteps;
for k = 1:numel(Sys)
  s = Sys(k);
  ok = ok && s.S==1/2 && s.charge==0 && isequal(s.Elements,{'O','H'});
  ok = ok && strcmp(s.Nucs,'O,H') && isequal(s.NucsIdx,[1 2]);
  ok = ok && isequal(size(s.g),[1 3]);
  for f = {'A','AFrame','Q','QFrame'}
    ok = ok && isequal(size(s.(f{1})),[2 3]);
  end
  ok = ok && abs(norm(s.xyz(2,:)-s.xyz(1,:))-rOH(k))<1e-4;
end
% Properties must change from step to step (each step is read separately)
for k = 2:numel(Sys)
  ok = ok && max(abs(Sys(k).A(2,:)-Sys(k-1).A(2,:)))>1;  % 1H, MHz
end
end

%-------------------------------------------------------------------------------
% Compare two spin system arrays step by step
function ok = comparescans(Sys1,Sys2)
T = @(v,ang) erot(ang).'*diag(v)*erot(ang);
maxdiff = @(X,Y) max(abs(X(:)-Y(:)));
ok = numel(Sys1)==numel(Sys2);
for k = 1:numel(Sys1)
  s1 = Sys1(k); s2 = Sys2(k);
  ok = ok && maxdiff(s1.xyz,s2.xyz)<1e-5;
  ok = ok && maxdiff(T(s1.g,s1.gFrame),T(s2.g,s2.gFrame))<1e-6;
  for n = 1:numel(s1.NucsIdx)
    ok = ok && maxdiff(T(s1.A(n,:),s1.AFrame(n,:)),T(s2.A(n,:),s2.AFrame(n,:)))<0.05;  % MHz
    ok = ok && maxdiff(T(s1.Q(n,:),s1.QFrame(n,:)),T(s2.Q(n,:),s2.QFrame(n,:)))<1e-3;  % MHz
  end
end
end
