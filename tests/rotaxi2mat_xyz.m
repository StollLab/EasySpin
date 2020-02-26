function ok = test()

% Explicitely check rotations around specific directions

R1 = rotaxi2mat('x',pi/2);
R2 = rotaxi2mat([1;0;0],pi/2);
Rref = [1 0 0; 0 0 1; 0 -1 0];
ok(1) = areequal(R1,Rref) && areequal(R2,Rref);

R1 = rotaxi2mat('y',pi/2);
R2 = rotaxi2mat([0;1;0],pi/2);
Rref = [0 0 -1; 0 1 0; 1 0 0];
ok(2) = areequal(R1,Rref) && areequal(R2,Rref);

R1 = rotaxi2mat('z',pi/2);
R2 = rotaxi2mat([0;0;1],pi/2);
Rref = [0 1 0; -1 0 0; 0 0 1];
ok(3) = areequal(R1,Rref) && areequal(R2,Rref);
