function ok = test()

B = rand(1,3)*400;
Sys = struct('S',1/2,'g',[2 3 4],'Nucs','63Cu','A',[50 50 350]);

H = ham(Sys,rand(3,1),'sparse');

[H0,Gx,Gy,Gz] = ham(Sys,[],'sparse');

ok = issparse(H) && issparse(H0) && issparse(Gx) && issparse(Gy) && issparse(Gz);
