function ok = test()

S = 1/2;
L = 1;
soc = 4545;

Sys.S = S;
Sys.L = L;
Sys.soc = soc;
H = soint(Sys);

[Sx,Sy,Sz] = sop([S L],'x1','y1','z1');
[Lx,Ly,Lz] = sop([S L],'x2','y2','z2');
Href = soc*(Sx*Lx + Sy*Ly + Sz*Lz);

ok = areequal(H,Href,1e-10,'abs');
