function ok = test()

% Make sure both open- and closed-phi grids work for C1 and Ci

n = 10;

sym = 'C1';
go = sphgrid(sym,n,'');   % open-phi grid
gc = sphgrid(sym,n,'c');  % closed-phi grid

ok(1) = gc.n==go.n+2*n-2;

sym = 'Ci';
go = sphgrid(sym,n,'');   % open-phi grid
gc = sphgrid(sym,n,'c');  % closed-phi grid

ok(2) = gc.n==go.n+n-1;
