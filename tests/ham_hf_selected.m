function ok = test()

rng(4545);

A{1,1} = rand(1,3);
A{1,2} = rand(1,3);
A{1,3} = rand(1,3);
A{2,1} = rand(1,3);
A{2,2} = rand(1,3);
A{2,3} = rand(1,3);

Sys.S = [1/2 1 1/2];
Sys.Nucs = '1H,1H';
Sys.A = [A{1,1} A{1,2} A{1,3}; A{2,1} A{2,2} A{2,3}];
Sys.ee = [1 1 1];

e = 2;
n = 2;

H1 = ham_hf(Sys,e,n);

spins = spinvec(Sys);
nEl = numel(Sys.S);

H2 = 0;
for c = 1:3
  H2 = H2 + A{e,n}(c)*sop(spins,[e c; nEl+n c]);
end

ok = areequal(H1,H2,1e-10,'abs');
