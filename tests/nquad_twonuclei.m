function ok = test()

% 2 nuclei

Sys = struct('S',1/2,'g',[2 2 2],'Nucs','2H,2H','A',[1 2 3; 4 5 6]);
Sys.Q = [-1 0 2; -1 0 2];
HQ = nquad(Sys,2);
QQ = [1.5 0 -0.5; 0 -1 0; -0.5 0 1.5];
QQ = repmat({QQ},1,6);
QQ = blkdiag(QQ{:});

ok = areequal(HQ,QQ,1e-10,'abs');
