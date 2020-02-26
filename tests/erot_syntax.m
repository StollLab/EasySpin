function ok = test()

% Test two input syntax variants

a = [1 2 3]/1.2343;
R1 = erot(a);
R2 = erot(a(1),a(2),a(3));

ok = areequal(R1,R2);
