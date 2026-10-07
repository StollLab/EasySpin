function ok = test()

% Symmetric matrices given as [xx yy zz xy xz yz] are expanded to full matrices

v = [1 2 3 0.1 0.2 0.3];
M = [1 0.1 0.2; 0.1 2 0.3; 0.2 0.3 3];
v2 = [4 5 6 -0.4 -0.5 -0.6];
M2 = [4 -0.4 -0.5; -0.4 5 -0.6; -0.5 -0.6 6];

% g, one electron
Sys = struct('S',1/2,'g',2+v/100);
Q = runprivate('validatespinsys',Sys);
ok(1) = Q.fullg && areequal(Q.g,2+M/100,1e-12,'abs');

% g, two electrons
Sys = struct('S',[1/2 1/2],'g',2+[v; v2]/100,'J',10);
Q = runprivate('validatespinsys',Sys);
ok(2) = Q.fullg && areequal(Q.g,2+[M; M2]/100,1e-12,'abs');

% D
Sys = struct('S',1,'D',v*100);
Q = runprivate('validatespinsys',Sys);
ok(3) = Q.fullD && areequal(Q.D,M*100,1e-12,'abs');

% ee
Sys = struct('S',[1/2 1/2],'ee',v);
Q = runprivate('validatespinsys',Sys);
ok(4) = Q.fullee && areequal(Q.ee,M,1e-12,'abs');

% A, one electron, two nuclei
Sys = struct('S',1/2,'Nucs','1H,1H','A',[v; v2]);
Q = runprivate('validatespinsys',Sys);
ok(5) = Q.fullA && areequal(Q.A,[M; M2],1e-12,'abs');

% A, two electrons, two nuclei
Sys = struct('S',[1/2 1/2],'J',10,'Nucs','1H,1H','A',[v v2; v2 v]);
Q = runprivate('validatespinsys',Sys);
ok(6) = Q.fullA && areequal(Q.A,[M M2; M2 M],1e-12,'abs');

% Q
Sys = struct('S',1/2,'Nucs','14N','A',[1 1 1],'Q',v-2);
Q = runprivate('validatespinsys',Sys);
ok(7) = Q.fullQ && areequal(Q.Q,M-2,1e-12,'abs');

% sigma
Sys = struct('S',1/2,'Nucs','1H','A',[1 1 1],'sigma',1+v*1e-3);
Q = runprivate('validatespinsys',Sys);
ok(8) = Q.fullsigma && areequal(Q.sigma,1+M*1e-3,1e-12,'abs');

% nn
Sys = struct('S',1/2,'Nucs','1H,1H','A',[1 1 1; 2 2 2],'nn',v);
Q = runprivate('validatespinsys',Sys);
ok(9) = Q.fullnn && areequal(Q.nn,M,1e-12,'abs');

% wrong sizes are still rejected
Sys = struct('S',1/2,'g',2+v(1:5)/100);
[~,err] = runprivate('validatespinsys',Sys);
ok(10) = ~isempty(err);

Sys = struct('S',1/2,'g',2+v.'/100);
[~,err] = runprivate('validatespinsys',Sys);
ok(11) = ~isempty(err);

% spin-zero nuclei are removed together with their AFrame/QFrame rows
Sys = struct('S',1/2,'Nucs','1H,12C','A',[v; v2],'Q',[v; v2]);
Q = runprivate('validatespinsys',Sys);
ok(12) = Q.nNuclei==1 && areequal(Q.A,M,1e-12,'abs') && ...
  isequal(size(Q.AFrame),[1 3]) && isequal(size(Q.QFrame),[1 3]);
