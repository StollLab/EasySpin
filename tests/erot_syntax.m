function ok = test()

% Test two input syntax variants

angles = [5 -15 76]*pi/180;
R1 = erot(angles);
R2 = erot(angles(1),angles(2),angles(3));

ok = areequal(R1,R2);
