function [err,data] = test(opt,olddata)

% Test whether sop returns sparse matrices when asked

Spins = [1/2 1];
Comp = '+z';

A = sop(Spins,Comp);
B = sop(Spins,Comp,'sparse');

ok = ~issparse(A) && issparse(B);
err = ~ok;

data = [];
