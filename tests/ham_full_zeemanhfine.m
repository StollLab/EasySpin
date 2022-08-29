function ok = test()

% ham vs zeeman/ham_hf

Sys = struct('S',1/2,'g',[2 3 4],'Nucs','63Cu','A',[50 50 350]);

B = rand(1,3)*400;
[F1,G1] = ham(Sys,B);

F = ham_nq(Sys) + ham_hf(Sys);
[eGx,eGy,eGz] = ham_ez(Sys);
[nGx,nGy,nGz] = ham_nz(Sys);
G = B(1)*(eGx+nGx) + B(2)*(eGy+nGy) + B(3)*(eGz+nGz);
G = G/norm(B);

ok = areequal(F1+norm(B)*G1,F+norm(B)*G,1e-10,'rel');
