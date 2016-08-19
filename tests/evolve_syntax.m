function [err,data] = test(opt,olddata)

herm = @(A)A+A';
randherm = @(n) herm(complex(rand(n),rand(n)));

N = 4;
Sig = randherm(N);
H = randherm(N);
M = expm(-1i*randherm(N));
D = Sig;
dt = 0.01;
n = 50;

% 1D experiments
td = evolve(Sig,D,H,n,dt,[1]);
td = evolve(Sig,D,H,n,dt,[1 1],M);
td = evolve(Sig,D,H,n,dt,[1 -1],M);
td = evolve(Sig,D,H,n,dt,[1 -1 1 -1],{M,M,M});
td = evolve(Sig,D,H,n,dt,[1 1 -1 -1],{M,M,M});
td = evolve(Sig,D,H,n,dt,[1 -1 -1 1],{M,M,M});

% 2D experiments
td = evolve(Sig,D,H,n,dt,[1 2],M);
td = evolve(Sig,D,H,n,dt,[1 2],{M});
td = evolve(Sig,D,H,n,dt,[1 1 2],{M,M});
td = evolve(Sig,D,H,n,dt,[1 -1 2],{M,M});
td = evolve(Sig,D,H,n,dt,[1 2 1],{M,M});
td = evolve(Sig,D,H,n,dt,[1 2 2 1],{M,M,M});
td = evolve(Sig,D,H,n,dt,[1 2 -2 1],{M,M,M});
td = evolve(Sig,D,H,n,dt,[1 1 -1 -1 2],{M,M,M,M});
td = evolve(Sig,D,H,n,dt,[1 -1 -1 1 2],{M,M,M,M});

err = 0;

data = [];
