function ok = test()

% Simple case 1 nucleus

Sys = struct('S',1/2,'g',2,'Nucs','2H','A',[1 2 3]);
Sys.Q = [4 6 8];
H = nquad(Sys);

Href = zeros(6,6);
Href([1 15 22 36]) = 13;
Href([8 29]) = 10;
Href([3 13 24 34]) = -1;

ok = areequal(H,Href,1e-10,'abs');
