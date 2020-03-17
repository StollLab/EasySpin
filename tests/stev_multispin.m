function ok = test()

% Stevens op. for multi-spin system: Check dimension

Spins = [3/2 3 5/2 3/2];

iSpin = 2;
k = 4;
q = +3;
Op = stev(Spins,[k,q,iSpin],'sparse');

Op0 = stev(Spins(iSpin),[k,q]);
n = 2*Spins+1;
npre = prod(n(1:iSpin-1));
npost = prod(n(iSpin+1:end));
Op0 = kron(kron(eye(npre),Op0),eye(npost));

ok = areequal(Op,Op0,1e-10,'abs');

