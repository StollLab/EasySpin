function err = test()

% Test two input syntax versions
%======================================================
a = [1 2 3]/1.2343;
R1 = erot(a);
R2 = erot(a(1),a(2),a(3));

err = ~areequal(R1,R2);
