function [err,data] = test(opt,olddata)

% Assert that sop returns sparse matrices when asked for
%-------------------------------------------------------

Spins = [1/2 1];
Comp1 = '+z';
Comp2 = 'x-';

% Case 1: single operator
A = sop(Spins,Comp1);
B = sop(Spins,Comp1,'sparse');

ok(1) = ~issparse(A) && issparse(B) && areequal(A,B,1e-12,'abs');

% Case 2: multiple operators
[A1,A2] = sop(Spins,Comp1,Comp2);
[B1,B2] = sop(Spins,Comp1,Comp2,'sparse');

ok(2) = ~issparse(A1) && ~issparse(A2) && ...
         issparse(B1) && issparse(B2) && ...
         areequal(A1,B1,1e-12,'abs') && areequal(A2,B2,1e-12,'abs');

err = ~ok;

data = [];
