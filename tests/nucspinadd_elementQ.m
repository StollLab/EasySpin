function ok = test()

% Adding nuclei to spin systems with quadrupole tensors, where nuclei are
% given as elements (as from orca2easyspin), or where Q is given in
% compact form and must be converted using the nuclear spin.

% Element-only Nucs with Q (as from orca2easyspin)
Sys.S = 1/2;
Sys.Nucs = 'C,N,H';
Sys.NucsIdx = [1 2 3];
Sys.A = [1 2 3; 4 5 6; 7 8 9];
Sys.Q = [0 0 0; -1 -1 2; 0 0 0];
Sys.QFrame = zeros(3,3);
Sys1 = nucspinadd(Sys,'14N',[1 2 3]);
ok(1) = strcmp(Sys1.Nucs,'C,N,H,14N') && isequal(Sys1.Q,[Sys.Q; 0 0 0]) && ...
  isequaln(Sys1.NucsIdx,[1 2 3 NaN]);

% Scalar Q for added element uses spin of Q reference isotope (63Cu, I=3/2)
Sys2 = nucspinadd(Sys,'Cu',[10 10 100],[],12);
ok(2) = areequal(Sys2.Q(end,:),[-1 -1 2],1e-12,'abs');

% Scalar Q for added nucleus uses its own spin, not those of prior nuclei
Sys = struct('S',1/2,'Nucs','63Cu,1H','A',[10 10 100; 1 2 3],'Q',[-1 -1 2; 0 0 0]);
Sys3 = nucspinadd(Sys,'14N',[1 2 3],[],4);  % I=1: 4/(4*1*1) = 1
ok(3) = areequal(Sys3.Q(end,:),[-1 -1 2],1e-12,'abs');

% Converting compact Q of existing nuclei: I=1/2 nucleus gives zero, not NaN
Sys = struct('S',1/2,'Nucs','1H,14N','A',[1 2 3; 4 5 6],'Q',[0; 4]);
Sys4 = nucspinadd(Sys,'2H',[1 1 1],[],[-0.1 -0.1 0.2]);
ok(4) = areequal(Sys4.Q,[0 0 0; -1 -1 2; -0.1 -0.1 0.2],1e-12,'abs');

% Spin system from ORCA property file
Sys = orca2easyspin('orcanew/v6.1.0/nitroxide.property.txt');
Sys5 = nucspinadd(Sys,'15N',[1 2 3]);
ok(5) = size(Sys5.Q,1)==numel(Sys.NucsIdx)+1 && isequal(Sys5.Q(1:end-1,:),Sys.Q);
