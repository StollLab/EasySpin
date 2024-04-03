function ok = test()

% Total transition intensity should be independent of orientation of k
% vector relative to B0, no matter what the crystal orientation and the
% polarization type and angle are.

Sys.S = 1;
Sys.D = 30e3*rand*10*[1 1/3];
Exp.SampleFrame = rand(1,3)*pi;
Opt.Transitions = [1 2; 1 3; 2 3];

for n=10:-1:1
  Exp.mwMode = {pi/2 rand*pi};
  [p(n,:),i(n,:)] = resfreqs_matrix(Sys,Exp,Opt);
end
ok(1) = all(abs(sum(i')-sum(i(1,:)))<1e-10);

for n=10:-1:1
  Exp.mwMode = {pi/2 pi/2};
  [p2(n,:),i2(n,:)] = resfreqs_matrix(Sys,Exp,Opt);
end
ok(2) = all(abs(sum(i2')-sum(i2(1,:)))<1e-10);

for n=10:-1:1
  Exp.mwMode = {rand*pi 'unpolarized'};
  [p3(n,:),i3(n,:)] = resfreqs_matrix(Sys,Exp,Opt);
end
ok(3) = all(abs(sum(i3')-sum(i3(1,:)))<1e-10);

for n = 10:-1:1
  Exp.mwMode = {rand*pi 'circular-'};
  [p4(n,:),i4(n,:)] = resfreqs_matrix(Sys,Exp,Opt);
end
ok(4) = all(abs(sum(i4')-sum(i4(1,:)))<1e-10);

for n = 10:-1:1
  Exp.mwMode = {rand*pi 'circular+'};
  [p5(n,:),i5(n,:)] = resfreqs_matrix(Sys,Exp,Opt);
end
ok(5) = all(abs(sum(i5')-sum(i5(1,:)))<1e-10);
