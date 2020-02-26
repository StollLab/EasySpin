function ok = test()

R0{1} = [0 0 -1; 0 1 0; 1 0 0];
R1{1} = erot(0,pi/2,0);

R0{2} = [0 1 0; 0 0 1; 1 0 0];
R1{2} = erot(0,pi/2,pi/2);

R0{3} = eye(3);
R1{3} = erot(0,0,0);

R0{4} = diag([1 -1 -1]);
R1{4} = erot(pi,-pi,0);

R0{5} = [0 1 0; -1 0 0; 0 0 1];
R1{5} = erot(pi/2,0,0);

for k = 1:numel(R0)
  ok(k) = areequal(R0{k},R1{k},1e-10,'abs');
end
