function ok = test()

% NucsIdx (atom indices, e.g. from orca2easyspin) must be trimmed
% together with Nucs, A, AFrame, Q, QFrame

Sys.S = 1/2;
Sys.Nucs = '13C,14N,17O,1H';
Sys.NucsIdx = [2 4 5 9];
Sys.A = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
Sys.AFrame = rand(4,3);
Sys.Q = [0 0 0; -1 -1 2; -2 -2 4; 0 0 0];
Sys.QFrame = rand(4,3);

Sys1 = nucspinkeep(Sys,2);
Sys2 = nucspinkeep(Sys,[1 4]);
Sys3 = nucspinrmv(Sys,[1 2 3 4]);

ok(1) = strcmp(Sys1.Nucs,'14N') && isequal(Sys1.NucsIdx,4) && ...
  isequal(Sys1.A,Sys.A(2,:)) && isequal(Sys1.Q,Sys.Q(2,:));
ok(2) = strcmp(Sys2.Nucs,'13C,1H') && isequal(Sys2.NucsIdx,[2 9]) && ...
  isequal(Sys2.AFrame,Sys.AFrame([1 4],:)) && isequal(Sys2.QFrame,Sys.QFrame([1 4],:));
ok(3) = ~isfield(Sys3,'Nucs') && ~isfield(Sys3,'NucsIdx');

% nucspinadd appends NaN, since the added nucleus has no atom index
Sys4 = nucspinadd(Sys2,'14N',[10 10 50]);
ok(4) = isequaln(Sys4.NucsIdx,[2 9 NaN]);
Sys5 = nucspinkeep(Sys4,[1 3]);
ok(5) = strcmp(Sys5.Nucs,'13C,14N') && isequaln(Sys5.NucsIdx,[2 NaN]);
Sys6 = nucspinadd(struct('S',1/2,'Nucs','1H','A',[1 2 3]),'14N',[10 10 50]);
ok(6) = ~isfield(Sys6,'NucsIdx');
