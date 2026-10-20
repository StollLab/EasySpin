function ok = test()

% Check that explicit spin indices are honored, for transition-selective
% operators, alpha/beta operators, and spin indices with more than one digit

tol = 1e-12;

% Transition-selective operators with explicit spin indices
A = sop([1/2 1/2],'x(1|2)2');
B = kron(eye(2),sop(1/2,'x(1|2)'));
ok(1) = areequal(A,B,tol,'abs');

A = sop([1/2 1 1/2],'x(1|3)2,z3');
B = kron(kron(eye(2),sop(1,'x(1|3)')),sop(1/2,'z'));
ok(2) = areequal(A,B,tol,'abs');

A = sop([3/2 1 1],'y(1|4)1,x(2|3)3');
B = kron(kron(sop(3/2,'y(1|4)'),eye(3)),sop(1,'x(2|3)'));
ok(3) = areequal(A,B,tol,'abs');

A = sop([1/2 1],'z(1|2)2,x1');
B = kron(sop(1/2,'x'),sop(1,'z(1|2)'));
ok(4) = areequal(A,B,tol,'abs');

% Alpha/beta operators with explicit spin indices
A = sop([1 1/2],'a2');
B = kron(eye(3),[1 0; 0 0]);
ok(5) = areequal(A,B,tol,'abs');

try
  sop([1/2 1],'b2');
  ok(6) = false;
catch
  ok(6) = true;
end

% Spin indices with more than one digit (spin-0 placeholders keep matrices small)
Spins = [1/2 zeros(1,10) 1];

A = sop(Spins,'x12,z1');
B = kron(sop(1/2,'z'),sop(1,'x'));
ok(7) = areequal(A,B,tol,'abs');

A = sop(Spins,'x(1|2)12,z1');
B = kron(sop(1/2,'z'),sop(1,'x(1|2)'));
ok(8) = areequal(A,B,tol,'abs');
