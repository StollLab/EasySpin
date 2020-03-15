function ok = test()

% Test multiple spins, with full kq specification

S = [2 3/2 1];
kq = [2 -1; 1 1; 2 0];
  
Op1 = isto(S(1),kq(1,:),'sparse');
Op2 = isto(S(2),kq(2,:),'sparse');
Op3 = isto(S(3),kq(3,:),'sparse');

Op = kron(Op1,Op2);
Op = kron(Op,Op3);

Op0 = isto(S,kq,'sparse');

ok = areequal(Op,Op0,1e-10,'abs');
