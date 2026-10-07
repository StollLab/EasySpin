function ok = test()

% Adding a nucleus with full matrices to a system with tilted principal-value
% tensors must give the same Hamiltonian as adding it with principal values

AFrame = [0.3 0.7 -0.4];
QFrame = [-0.2 0.5 1.1];
Sys = struct('S',1/2,'Nucs','14N','A',[10 20 30],'AFrame',AFrame,'Q',[-1 -1 2],'QFrame',QFrame);

Apv = [1 2 4];
Qpv = [-0.2 -0.1 0.3];

% new nucleus with principal values
Sys1 = nucspinadd(Sys,'2H',Apv,[0 0 0],Qpv,[0 0 0]);

% new nucleus with full matrices
Sys2 = nucspinadd(Sys,'2H',diag(Apv),[],diag(Qpv),[]);

H1 = ham(Sys1,[100 200 300]);
H2 = ham(Sys2,[100 200 300]);
ok(1) = areequal(H1,H2,1e-10,'rel');

% adding a tilted principal-value nucleus to a system with full matrices
SysF = struct('S',1/2,'Nucs','14N','A',diag([10 20 30]),'Q',diag([-1 -1 2]));
Sys3 = nucspinadd(SysF,'2H',Apv,AFrame,Qpv,QFrame);
RA = erot(AFrame);
RQ = erot(QFrame);
Sys4 = struct('S',1/2,'Nucs','14N,2H',...
  'A',[SysF.A; RA.'*diag(Apv)*RA],'Q',[SysF.Q; RQ.'*diag(Qpv)*RQ]);
ok(end+1) = areequal(ham(Sys3,[100 200 300]),ham(Sys4,[100 200 300]),1e-10,'rel');

% full matrix with nonzero frame is rejected
try
  nucspinadd(Sys,'2H',diag(Apv),[0 1 0]);
  ok(end+1) = false;
catch
  ok(end+1) = true;
end
